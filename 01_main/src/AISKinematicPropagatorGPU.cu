// AISKinematicPropagatorGPU.cu
//
// GPU acceleration for kinematic (free-flight) propagation in ultrafast mode.
//
// Design principles:
//   - Compiled only when -DUSE_CUDA=ON.
//   - CPU paths in AISKinematicPropagator.cc are completely unchanged.
//   - Supports zero_pot and linear_pot (both have constant acceleration
//     so the trajectory update is an exact analytic expression, not ODE).
//   - Falls back silently to CPU for quadratic_pot (position-dependent force).
//   - The __float128 quad-phase update (phaseQuad -= omegaSr87 * dt for
//     excited-state wavepackets) stays on the CPU — it's a uniform scalar
//     subtraction done in a tight loop, negligible compared to the ODE cost.

#include "AISKinematicPropagatorGPU.hh"
#include "AISKinematicPropagator.hh"
#include "AISAtomEnsemble.hh"
#include "AISWavePacket.hh"
#include "AISAtom.hh"
#include "AISConstants.hh"

#include <cuda_runtime.h>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <memory>

// ---------------------------------------------------------------------------
// Error helper (same pattern as pulse propagator GPU file)
// ---------------------------------------------------------------------------
static void cudaCheck(cudaError_t err, const char* file, int line)
{
    if (err != cudaSuccess) {
        std::ostringstream oss;
        oss << "CUDA error at " << file << ":" << line
            << " — " << cudaGetErrorString(err);
        throw std::runtime_error(oss.str());
    }
}
#define CUDA_CHECK(x) cudaCheck((x), __FILE__, __LINE__)

// ---------------------------------------------------------------------------
// Kernel: exact analytic trajectory update for constant acceleration.
//
// For zero_pot:    a = {0, 0, 0}      → pos += vel*dt, vel unchanged
// For linear_pot:  a = {0, 0, -g}     → uniform gravity
//
// In-place update of pos and vel arrays.
// ---------------------------------------------------------------------------
__global__ void kinematicPropagateKernel(
    double* __restrict__ pos_x,
    double* __restrict__ pos_y,
    double* __restrict__ pos_z,
    double* __restrict__ vel_x,
    double* __restrict__ vel_y,
    double* __restrict__ vel_z,
    double ax, double ay, double az,
    double dt,
    int N)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    double vx = vel_x[i], vy = vel_y[i], vz = vel_z[i];
    double half_dt2 = 0.5 * dt * dt;

    pos_x[i] += vx * dt + ax * half_dt2;
    pos_y[i] += vy * dt + ay * half_dt2;
    pos_z[i] += vz * dt + az * half_dt2;

    vel_x[i] += ax * dt;
    vel_y[i] += ay * dt;
    vel_z[i] += az * dt;
}

// ---------------------------------------------------------------------------
// launchKinematicKernel — host-side wrapper (declared in the .hh)
// ---------------------------------------------------------------------------
void launchKinematicKernel(
    double* pos_x, double* pos_y, double* pos_z,
    double* vel_x, double* vel_y, double* vel_z,
    double ax, double ay, double az,
    double dt,
    int N)
{
    if (N == 0) return;

    size_t sz = N * sizeof(double);

    double *d_px, *d_py, *d_pz, *d_vx, *d_vy, *d_vz;
    CUDA_CHECK(cudaMalloc(&d_px, sz)); CUDA_CHECK(cudaMalloc(&d_py, sz));
    CUDA_CHECK(cudaMalloc(&d_pz, sz)); CUDA_CHECK(cudaMalloc(&d_vx, sz));
    CUDA_CHECK(cudaMalloc(&d_vy, sz)); CUDA_CHECK(cudaMalloc(&d_vz, sz));

    // H→D
    CUDA_CHECK(cudaMemcpy(d_px, pos_x, sz, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_py, pos_y, sz, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_pz, pos_z, sz, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_vx, vel_x, sz, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_vy, vel_y, sz, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_vz, vel_z, sz, cudaMemcpyHostToDevice));

    int blockSize = 256;
    int gridSize  = (N + blockSize - 1) / blockSize;
    kinematicPropagateKernel<<<gridSize, blockSize>>>(
        d_px, d_py, d_pz, d_vx, d_vy, d_vz, ax, ay, az, dt, N);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    // D→H
    CUDA_CHECK(cudaMemcpy(pos_x, d_px, sz, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(pos_y, d_py, sz, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(pos_z, d_pz, sz, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(vel_x, d_vx, sz, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(vel_y, d_vy, sz, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(vel_z, d_vz, sz, cudaMemcpyDeviceToHost));

    cudaFree(d_px); cudaFree(d_py); cudaFree(d_pz);
    cudaFree(d_vx); cudaFree(d_vy); cudaFree(d_vz);
}

// ---------------------------------------------------------------------------
// AISKinematicPropagator::PropagateEnsembleGPU
//
// Replaces PropagateEnsemble when USE_CUDA is on and potentialType is
// zero_pot or linear_pot (constant acceleration — analytic trajectory).
//
// Structure:
//   1. Flatten pos/vel of all wavepackets across all atoms.
//   2. GPU kernel: analytic pos/vel update.
//   3. Write pos/vel back to wavepackets.
//   4. CPU: quad-phase update for state==1 wavepackets (uniform scalar op).
//   5. CPU: update wavepacket time stamp.
// ---------------------------------------------------------------------------
void AISKinematicPropagator::PropagateEnsembleGPU(
    std::unique_ptr<AISAtomEnsemble>& atomEnsemble, __float128 t1)
{
    int nAtoms = atomEnsemble->GetNumberOfAtoms();

    // ----------------------------------------------------------------
    // Step 0: constants for this propagation step
    // ----------------------------------------------------------------
    // t0 is read per-wavepacket (they all share the same t0 after a pulse,
    // but we use the first atom's first wavepacket as reference — all are equal).
    __float128 t0_128 = atomEnsemble->GetAtom(0)->GetWavePacket(0)->GetTime();
    double t0d = static_cast<double>(t0_128);
    double t1d = static_cast<double>(t1);
    double dt  = t1d - t0d;

    // Constant acceleration: -dUdx / massSr87
    //   zero_pot:   a = {0, 0, 0}
    //   linear_pot: dUdx = {0, 0, massSr87*g}  →  a = {0, 0, -g}
    double ax = constAcceleration[0];
    double ay = constAcceleration[1];
    double az = constAcceleration[2];

    // quad-phase shift for state==1 wavepackets: phaseQuad -= omegaSr87 * dt
    // Computed once on CPU; same for every excited-state wavepacket.
    __float128 phaseQuadDelta = omegaSr87 * (t1 - t0_128);

    // ----------------------------------------------------------------
    // Step 1: flatten all wavepackets into host arrays
    // ----------------------------------------------------------------
    // Count total wavepackets
    int totalWP = 0;
    for (int a = 0; a < nAtoms; ++a)
        totalWP += atomEnsemble->GetAtom(a)->GetNumberOfWavePackets();

    std::vector<double> h_px(totalWP), h_py(totalWP), h_pz(totalWP);
    std::vector<double> h_vx(totalWP), h_vy(totalWP), h_vz(totalWP);

    {
        int idx = 0;
        for (int a = 0; a < nAtoms; ++a) {
            auto& atom = atomEnsemble->GetAtom(a);
            int nWP = atom->GetNumberOfWavePackets();
            for (int j = 0; j < nWP; ++j, ++idx) {
                auto& wp = atom->GetWavePacket(j);
                auto pos = wp->GetPosition();
                auto vel = wp->GetVelocity();
                h_px[idx] = pos[0]; h_py[idx] = pos[1]; h_pz[idx] = pos[2];
                h_vx[idx] = vel[0]; h_vy[idx] = vel[1]; h_vz[idx] = vel[2];
            }
        }
    }

    // ----------------------------------------------------------------
    // Step 2: GPU kernel — analytic trajectory update
    // ----------------------------------------------------------------
    launchKinematicKernel(
        h_px.data(), h_py.data(), h_pz.data(),
        h_vx.data(), h_vy.data(), h_vz.data(),
        ax, ay, az, dt, totalWP);

    // ----------------------------------------------------------------
    // Step 3: write updated pos/vel back + quad phase update (CPU) + timestamp
    // ----------------------------------------------------------------
    {
        int idx = 0;
        for (int a = 0; a < nAtoms; ++a) {
            auto& atom = atomEnsemble->GetAtom(a);
            int nWP = atom->GetNumberOfWavePackets();
            for (int j = 0; j < nWP; ++j, ++idx) {
                auto& wp = atom->GetWavePacket(j);

                wp->SetPosition({h_px[idx], h_py[idx], h_pz[idx]});
                wp->SetVelocity({h_vx[idx], h_vy[idx], h_vz[idx]});

                // phaseDouble is unchanged in ultrafast mode (no action phase)
                // phaseQuad: same scalar subtracted from all state==1 wavepackets
                if (wp->GetState() == 1)
                    wp->SetPhaseQuad(wp->GetPhaseQuad() - phaseQuadDelta);

                wp->SetTime(t1);
            }
        }
    }
}
