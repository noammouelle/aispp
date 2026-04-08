// AISPulsePropagatorGPU.cu
//
// CUDA implementation of the UltraFast U3 pulse propagator.
//
// Design principles:
//   - This file is compiled ONLY when -DUSE_CUDA=ON is passed to CMake.
//   - The CPU code paths in AISPulsePropagator.cc are completely unchanged.
//   - The only entry point visible to the rest of the codebase is
//     AISPulsePropagator::PropagateEnsembleGPU(), which is called from
//     PropagateEnsemble() under a #ifdef USE_CUDA guard.
//   - Zero __float128 arithmetic appears inside the kernel.  The one quad
//     subtraction (omega - omegaSr87) is performed on the CPU before the
//     kernel launch and the double result is passed in GPUBeamParams.
//
// Parallelism strategy:
//   One CUDA thread per wavepacket.  Because atoms are independent and
//   wavepackets within an atom are independent at the U3 stage, this is
//   embarrassingly parallel with no inter-thread communication.

#include "AISPulsePropagatorGPU.hh"
#include "AISPulsePropagator.hh"
#include "AISAtomEnsemble.hh"
#include "AISWavePacket.hh"
#include "AISAtom.hh"
#include "AISLaserBeam.hh"
#include "AISConstants.hh"

#include <cuda_runtime.h>
#include <omp.h>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <memory>
#include <cmath>
#include <string>

// ---------------------------------------------------------------------------
// Helpers
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
// Physical constants — replicated as constexpr for device code.
// Values match AISConstants.hh exactly (computed at double precision).
// ---------------------------------------------------------------------------

static constexpr double GPU_PI          = 3.141592653589793;
static constexpr double GPU_HBAR        = 1.054571817646157e-34; // J·s
static constexpr double GPU_MASS_SR87   = 1.443155904223061e-25; // kg
static constexpr double GPU_LAMBDA_SR87 = 6.984457096127544e-07; // m  (c/f_Sr87)

// ---------------------------------------------------------------------------
// Device helper: Gaussian intensity envelope
//   pos is ALREADY shifted by focal length (caller's responsibility).
// ---------------------------------------------------------------------------
__device__ static double gaussianEnvelopeGPU(double rx, double ry, double rz, double w0)
{
    double r2 = rx * rx + ry * ry;
    double zR = GPU_PI * w0 * w0 / GPU_LAMBDA_SR87;
    double wz2 = w0 * w0 * (1.0 + (rz / zR) * (rz / zR));
    return exp(-r2 / wz2) * (w0 / sqrt(wz2));
}

// ---------------------------------------------------------------------------
// Device helper: effective Rabi frequency at atom position
// ---------------------------------------------------------------------------
__device__ static double getRabiFreqGPU(double px, double py, double pz,
                                        const GPUBeamParams bp)
{
    if (bp.beamTypeFlat) return bp.rabiFreq;

    // Gaussian beam: shift z by focal length (sign depends on beam direction)
    double zShifted = bp.focalLength + bp.kzSign * pz;
    return bp.rabiFreq * gaussianEnvelopeGPU(px, py, zShifted, bp.w0);
}

// ---------------------------------------------------------------------------
// Kernel: one thread per (atom, wavepacket) pair.
//
// Implements ApplyU3UltraFast exactly as in AISPulsePropagator.cc but for
// all wavepackets simultaneously.
//
// Output convention:
//   *_out[i] is the result for input wavepacket i.
//   wavepacket-0 ("keep-state"):  amp0_out, dphase0_out
//   wavepacket-1 ("flip-state"):  amp1_out, dphase1_out
//   dphase = arg(S_matrix_element)  — added to phaseDouble on the CPU side.
// ---------------------------------------------------------------------------
__global__ void applyU3UltraFastKernel(
    const double* __restrict__ pos_x,
    const double* __restrict__ pos_y,
    const double* __restrict__ pos_z,
    const double* __restrict__ vel_x,
    const double* __restrict__ vel_y,
    const double* __restrict__ vel_z,
    const double* __restrict__ amp_in,
    const int*    __restrict__ state_in,
    double* __restrict__ amp0_out,
    double* __restrict__ dphase0_out,
    double* __restrict__ amp1_out,
    double* __restrict__ dphase1_out,
    GPUBeamParams bp,
    int N)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    double px = pos_x[i], py = pos_y[i], pz = pos_z[i];
    double vx = vel_x[i], vy = vel_y[i], vz = vel_z[i];
    int    st = state_in[i];
    double amp = amp_in[i];

    // --- Detuning (mirrors getDeltaUltraFast, gradPhi==0 so kPrime==kT1) ---
    // Doppler term: dot(kT1, vel)
    double dopplerTerm = bp.kT1_x * vx + bp.kT1_y * vy + bp.kT1_z * vz;
    // Wavevector chirp term: 0.5 * dot(kChirp, pos)
    double kChirpTerm  = 0.5 * (bp.kChirp_x * px + bp.kChirp_y * py + bp.kChirp_z * pz);
    // bp.recoilTerm and bp.dOmegaScalar and bp.omegaChirpTerm are global constants
    double detuning = dopplerTerm - bp.dOmegaScalar + bp.recoilTerm
                      + kChirpTerm + bp.omegaChirpTerm;

    // --- Effective Rabi frequency ---
    double Omega = getRabiFreqGPU(px, py, pz, bp);

    // --- 2-level SU(2) unitary (analytic, same as ApplyU3UltraFast) ---
    double OmegaEff = sqrt(Omega * Omega + detuning * detuning);
    double halfDt   = 0.5 * (bp.t1 - bp.t0);
    double trigArg  = halfDt * OmegaEff;
    double expPhase = -halfDt * detuning;  // phase of expFactor = exp(i*expPhase)

    double cosT = cos(trigArg);
    double sinT = sin(trigArg);
    double cosE = cos(expPhase);
    double sinE = sin(expPhase);

    // Guard against zero OmegaEff (degenerate case: no drive, no detuning)
    double sinRatio = (OmegaEff > 1e-30) ? (Omega    / OmegaEff) * sinT : 0.0;
    double detRatio = (OmegaEff > 1e-30) ? (detuning / OmegaEff) * sinT : 0.0;

    // Matrix elements — keep-state (Skk) and flip-state (Sflip)
    // state==0: Skk = Sgg, Sflip = Seg
    // state==1: Skk = See, Sflip = Sge
    // Note: Seg == Sge (same formula), only Skk differs between states.
    double Skk_re, Skk_im;
    if (st == 0) {
        // Sgg = expFactor * (cosT + i*detRatio)
        Skk_re = cosE * cosT - sinE * detRatio;
        Skk_im = sinE * cosT + cosE * detRatio;
    } else {
        // See = expFactor * (cosT - i*detRatio)
        Skk_re = cosE * cosT + sinE * detRatio;
        Skk_im = sinE * cosT - cosE * detRatio;
    }
    // Seg = Sge = expFactor * (-i * sinRatio)
    double Sflip_re =  sinE * sinRatio;
    double Sflip_im = -cosE * sinRatio;

    amp0_out[i]    = amp * hypot(Skk_re,   Skk_im);
    dphase0_out[i] = atan2(Skk_im,   Skk_re);
    amp1_out[i]    = amp * hypot(Sflip_re, Sflip_im);
    dphase1_out[i] = atan2(Sflip_im, Sflip_re);
}

// ---------------------------------------------------------------------------
// Persistent device-side storage for pulse propagator.
//
// All 11 double arrays and the int state array are packed into a single
// cudaMalloc.  The block is grown (never shrunk) as N increases, so
// malloc/free overhead is paid at most O(log N_max) times over the entire
// simulation rather than once per pulse step.
// ---------------------------------------------------------------------------
struct PulseDeviceBuffers {
    void*  block    = nullptr;
    int    capacity = 0;

    double *d_px, *d_py, *d_pz;
    double *d_vx, *d_vy, *d_vz;
    double *d_amp;
    double *d_amp0, *d_dph0, *d_amp1, *d_dph1;
    int    *d_state;

    void ensure(int N) {
        if (N <= capacity) return;
        if (block) { cudaFree(block); block = nullptr; }

        // 11 double arrays + state ints packed after the doubles (double-aligned).
        size_t totalBytes = 11 * (size_t)N * sizeof(double)
                          +      (size_t)N * sizeof(int);
        CUDA_CHECK(cudaMalloc(&block, totalBytes));

        double* base = static_cast<double*>(block);
        d_px    = base +  0 * N;
        d_py    = base +  1 * N;
        d_pz    = base +  2 * N;
        d_vx    = base +  3 * N;
        d_vy    = base +  4 * N;
        d_vz    = base +  5 * N;
        d_amp   = base +  6 * N;
        d_amp0  = base +  7 * N;
        d_dph0  = base +  8 * N;
        d_amp1  = base +  9 * N;
        d_dph1  = base + 10 * N;
        d_state = reinterpret_cast<int*>(base + 11 * N);  // double-aligned boundary
        capacity = N;
    }

    ~PulseDeviceBuffers() { if (block) cudaFree(block); }
};

static PulseDeviceBuffers s_pulseBufs;

// ---------------------------------------------------------------------------
// launchU3UltraFastKernel — host-side wrapper (declared in the .hh)
//
// Uses persistent device buffers to avoid per-call cudaMalloc/cudaFree.
// Input arrays are host pointers; output arrays are host pointers.
// ---------------------------------------------------------------------------
void launchU3UltraFastKernel(
    const double* pos_x, const double* pos_y, const double* pos_z,
    const double* vel_x, const double* vel_y, const double* vel_z,
    const double* amp_in,
    const int*    state_in,
    double* amp0_out,   double* dphase0_out,
    double* amp1_out,   double* dphase1_out,
    const GPUBeamParams& bp,
    int N)
{
    if (N == 0) return;

    // --- Ensure device buffers are large enough (no-op if N hasn't grown) ---
    s_pulseBufs.ensure(N);

    size_t szD = N * sizeof(double);
    size_t szI = N * sizeof(int);

    // --- Copy inputs H→D ---
    CUDA_CHECK(cudaMemcpy(s_pulseBufs.d_px,    pos_x,    szD, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(s_pulseBufs.d_py,    pos_y,    szD, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(s_pulseBufs.d_pz,    pos_z,    szD, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(s_pulseBufs.d_vx,    vel_x,    szD, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(s_pulseBufs.d_vy,    vel_y,    szD, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(s_pulseBufs.d_vz,    vel_z,    szD, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(s_pulseBufs.d_amp,   amp_in,   szD, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(s_pulseBufs.d_state, state_in, szI, cudaMemcpyHostToDevice));

    // --- Launch kernel ---
    int blockSize = 256;
    int gridSize  = (N + blockSize - 1) / blockSize;
    applyU3UltraFastKernel<<<gridSize, blockSize>>>(
        s_pulseBufs.d_px,  s_pulseBufs.d_py,  s_pulseBufs.d_pz,
        s_pulseBufs.d_vx,  s_pulseBufs.d_vy,  s_pulseBufs.d_vz,
        s_pulseBufs.d_amp, s_pulseBufs.d_state,
        s_pulseBufs.d_amp0, s_pulseBufs.d_dph0,
        s_pulseBufs.d_amp1, s_pulseBufs.d_dph1,
        bp, N);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    // --- Copy outputs D→H ---
    CUDA_CHECK(cudaMemcpy(amp0_out,    s_pulseBufs.d_amp0, szD, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(dphase0_out, s_pulseBufs.d_dph0, szD, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(amp1_out,    s_pulseBufs.d_amp1, szD, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(dphase1_out, s_pulseBufs.d_dph1, szD, cudaMemcpyDeviceToHost));
}

// ---------------------------------------------------------------------------
// AISPulsePropagator::PropagateEnsembleGPU
//
// Replaces PropagateEnsemble for ultrafast mode.  Structure mirrors
// PropagateAtom() exactly — only the U3UltraFast step runs on the GPU;
// U2, U2†, U1, path-selection, cutoff all remain on the CPU and are
// called via the same existing methods.
// ---------------------------------------------------------------------------
void AISPulsePropagator::PropagateEnsembleGPU(
    std::unique_ptr<AISAtomEnsemble>& atomEnsemble)
{
    int nAtoms = atomEnsemble->GetNumberOfAtoms();

    // ----------------------------------------------------------------
    // Step 0: Pre-compute global scalars on CPU.
    //         This is the ONLY __float128 arithmetic in this function.
    // ----------------------------------------------------------------
    double t0d = this->initTimeDouble;
    double t1d = this->finalTimeDouble;

    // k vector at t1 (chirped)
    doubleThreeVector kT1    = laserBeam->GetK(this->finalTime);
    doubleThreeVector kChirp = laserBeam->GetKChirp();

    // THE key computation: quad subtraction done once, result stored as double.
    // This is safe because |omega - omegaSr87| ~ MHz << omega ~ 10^15 rad/s,
    // so the result fits well within double precision range.
    __float128 dOmega128  = laserBeam->GetOmega(0.0q) - omegaSr87;
    double dOmegaScalar   = static_cast<double>(dOmega128);

    // Frequency chirp term: -omegaChirp * t1 (same for all atoms)
    double omegaChirpTerm = -static_cast<double>(laserBeam->GetFrequencyChirp()) * t1d;

    // Recoil term: |kT1|^2 * hbar/(2*m)
    // gradPhi == {0,0,0} currently, so kPrime == kT1 for all atoms.
    double recoilTerm = dotProduct(kT1, kT1) * hbar / (2.0 * massSr87);

    // Beam parameters struct for kernel
    GPUBeamParams bp;
    bp.kT1_x = kT1[0];    bp.kT1_y = kT1[1];    bp.kT1_z = kT1[2];
    bp.kChirp_x = kChirp[0]; bp.kChirp_y = kChirp[1]; bp.kChirp_z = kChirp[2];
    bp.rabiFreq       = laserBeam->GetW0() > 0 ? laserBeam->GetW0() : 0; // placeholder
    bp.w0             = laserBeam->GetW0();
    bp.focalLength    = laserBeam->GetFocalLength();
    bp.dOmegaScalar   = dOmegaScalar;
    bp.omegaChirpTerm = omegaChirpTerm;
    bp.recoilTerm     = recoilTerm;
    bp.t0             = t0d;
    bp.t1             = t1d;
    bp.kzSign         = (kT1[2] >= 0.0) ? 1 : -1;
    bp.beamTypeFlat   = (laserBeam->GetBeamType() == "flat_square") ? 1 : 0;

    // Retrieve the actual central Rabi frequency.
    // GetRabiFreq takes a position; use origin as reference — the envelope
    // factor is applied per-atom inside the kernel via getRabiFreqGPU.
    // We need the *central* rabiFreq field of the beam object.
    // Since it's private in AISLaserBeam we recover it by evaluating at
    // the beam axis (r=0, z=focalLength → envelope=1 for Gaussian, or 1 for flat).
    {
        doubleThreeVector axisPos = {0.0, 0.0, 0.0};
        // For gaussian the envelope at r=0, z=shifted-focal-length = 1 by construction
        // when shiftedZ == 0; but shiftedZ = focalLength + kzSign*0 = focalLength,
        // and w(focalLength) may not equal w0.  Instead, just query at the centre of
        // the beam: pos such that the shifted z = 0 (beam waist).
        // Simpler: call GetRabiFreq with a position that gives envelope=1.
        // For flat_square this is any position. For Gaussian this is the waist.
        // kzSign * pos[2] = -focalLength => pos[2] = -focalLength * kzSign
        doubleThreeVector waistPos = {0.0, 0.0, -(double)bp.kzSign * bp.focalLength};
        // at this position shiftedZ = focalLength + kzSign * (-kzSign*focalLength) = 0
        // => envelope = w0/w(0) = 1. So GetRabiFreq returns rabiFreq_central.
        bp.rabiFreq = laserBeam->GetRabiFreq(waistPos, this->initTime, this->finalTime);
    }

    // ----------------------------------------------------------------
    // Step 1: Apply U2 — parallel over atoms.
    //   U2 only modifies state==1 wavepackets and uses quad arithmetic.
    //   Also set posStar/velStar which are used later in U2Dagger/U1.
    //   Each atom is fully independent so this is embarrassingly parallel.
    // ----------------------------------------------------------------
    #pragma omp parallel for schedule(static)
    for (int a = 0; a < nAtoms; ++a) {
        auto& atom = atomEnsemble->GetAtom(a);
        for (int j = 0; j < atom->GetNumberOfWavePackets(); ++j) {
            auto& wp = atom->GetWavePacket(j);
            wp->SetPosStar(wp->GetPosition());
            wp->SetVelStar(wp->GetVelocity());
            ApplyU2(wp, this->initTime);
        }
    }

    // ----------------------------------------------------------------
    // Step 2: Compute per-atom WP offsets, then flatten into SoA arrays.
    //   Offset array lets the flatten and unpack loops run in parallel
    //   without a shared sequential counter.
    //   WP counts are gathered in parallel (OMP) to avoid serial pointer-
    //   chasing over the atom list; the prefix sum over plain ints is fast.
    // ----------------------------------------------------------------
    std::vector<int> wpCounts(nAtoms);
    #pragma omp parallel for schedule(static)
    for (int a = 0; a < nAtoms; ++a)
        wpCounts[a] = atomEnsemble->GetAtom(a)->GetNumberOfWavePackets();

    std::vector<int> atomOffset(nAtoms + 1, 0);
    for (int a = 0; a < nAtoms; ++a)
        atomOffset[a + 1] = atomOffset[a] + wpCounts[a];
    int N = atomOffset[nAtoms];

    std::vector<double> h_px(N), h_py(N), h_pz(N);
    std::vector<double> h_vx(N), h_vy(N), h_vz(N);
    std::vector<double> h_amp(N);
    std::vector<int>    h_state(N);

    #pragma omp parallel for schedule(static)
    for (int a = 0; a < nAtoms; ++a) {
        int base = atomOffset[a];
        auto& atom = atomEnsemble->GetAtom(a);
        for (int j = 0; j < atom->GetNumberOfWavePackets(); ++j) {
            int i = base + j;
            auto& wp = atom->GetWavePacket(j);
            auto pos = wp->GetPosition();
            auto vel = wp->GetVelocity();
            h_px[i] = pos[0]; h_py[i] = pos[1]; h_pz[i] = pos[2];
            h_vx[i] = vel[0]; h_vy[i] = vel[1]; h_vz[i] = vel[2];
            h_amp[i]   = wp->GetAmplitude();
            h_state[i] = wp->GetState();
        }
    }

    // ----------------------------------------------------------------
    // Step 3: GPU kernel — computes U3UltraFast for all wavepackets.
    // ----------------------------------------------------------------
    std::vector<double> h_amp0(N), h_dph0(N), h_amp1(N), h_dph1(N);

    launchU3UltraFastKernel(
        h_px.data(), h_py.data(), h_pz.data(),
        h_vx.data(), h_vy.data(), h_vz.data(),
        h_amp.data(), h_state.data(),
        h_amp0.data(), h_dph0.data(),
        h_amp1.data(), h_dph1.data(),
        bp, N);

    // ----------------------------------------------------------------
    // Step 4+5: Unpack GPU results, create daughter wavepackets, apply
    //   post-processing and inverse transforms — parallel over atoms.
    //   Each atom's newWavePackets vector is local, so no data races.
    //
    // Fast-path optimisation: when no MC branching, path selection, or
    // det-vol selection is active (the common case), we skip the per-WP
    // tempVec allocation entirely and inline the amplitude cutoff.
    // A thread_local accumulator retains its heap buffer across atoms
    // handled by the same thread, eliminating repeated resize costs.
    // ----------------------------------------------------------------
    const bool fastPath = !this->useMcBranching &&
                          !this->usePathSelection &&
                          !this->useDetVolSelection;

    #pragma omp parallel for schedule(dynamic, 64)
    for (int a = 0; a < nAtoms; ++a) {
        int base   = atomOffset[a];
        auto& atom = atomEnsemble->GetAtom(a);
        int nWP    = atom->GetNumberOfWavePackets();

        // Thread-local accumulator: clear() keeps the heap buffer alive
        // across atoms on the same thread, amortising resize allocations.
        thread_local wavePacketVector tl_newWPs;
        tl_newWPs.clear();

        for (int j = 0; j < nWP; ++j) {
            int flatIdx = base + j;
            auto& wp0 = atom->GetWavePacket(j);

            // Build wavepacket1 (the "flip-state" daughter)
            std::unique_ptr<AISWavePacket> wp1(new AISWavePacket());
            wp1->SetPos0(wp0->GetPos0());
            wp1->SetVel0(wp0->GetVel0());
            wp1->SetT0(wp0->GetT0());

            // Bookkeeping copied from wavepacket0 before it is modified
            double      curPhaseD = wp0->GetPhaseDouble();
            __float128  curPhaseQ = wp0->GetPhaseQuad();
            auto pos              = wp0->GetPosition();
            auto vel              = wp0->GetVelocity();
            auto posStar          = wp0->GetPosStar();
            auto velStar          = wp0->GetVelStar();
            std::string path      = wp0->GetPath();
            auto detectPaths      = wp0->GetDetectablePaths();
            bool willInterf       = wp0->GetWillInterfere();
            int  st               = wp0->GetState();

            char keepSuffix = (st == 0) ? '0' : '1';
            char flipSuffix = (st == 0) ? '1' : '0';

            // --- Update wavepacket0 (keep-state) ---
            wp0->SetAmplitude(h_amp0[flatIdx]);
            wp0->SetPhaseDouble(curPhaseD + h_dph0[flatIdx]);
            wp0->SetPath(path + keepSuffix);
            wp0->SetPosStar(posStar);
            wp0->SetVelStar(velStar);
            wp0->SetDetectablePaths(detectPaths);
            wp0->SetWillInterfere(willInterf);

            // --- Populate wavepacket1 (flip-state) ---
            wp1->SetAmplitude(h_amp1[flatIdx]);
            wp1->SetPhaseDouble(curPhaseD + h_dph1[flatIdx]);
            wp1->SetPhaseQuad(curPhaseQ);
            wp1->SetState(1 - st);
            wp1->SetPosition(pos);
            wp1->SetVelocity(vel);
            wp1->SetPath(path + flipSuffix);
            wp1->SetPosStar(posStar);
            wp1->SetVelStar(velStar);
            wp1->SetDetectablePaths(detectPaths);
            wp1->SetWillInterfere(willInterf);

            if (fastPath) {
                // Inline amplitude cutoff — avoids two 'new wavePacketVector'
                // per WP (one for tempVec, one inside ApplyCutoff).
                if (std::abs(h_amp0[flatIdx]) > this->amplitudeThreshold)
                    tl_newWPs.push_back(std::move(wp0));
                if (std::abs(h_amp1[flatIdx]) > this->amplitudeThreshold)
                    tl_newWPs.push_back(std::move(wp1));
            } else {
                // Slow path: Apply* filters require unique_ptr<wavePacketVector>.
                std::unique_ptr<wavePacketVector> tempVec(new wavePacketVector);
                tempVec->push_back(std::move(wp0));
                tempVec->push_back(std::move(wp1));
                if (this->useMcBranching)     ApplyMCBranching(tempVec);
                if (this->usePathSelection)   ApplyPathSelection(tempVec, pathsToSimulate);
                if (this->useDetVolSelection) ApplyDetVolSelection(tempVec);
                ApplyCutoff(tempVec);
                for (auto& wp : *tempVec)
                    tl_newWPs.push_back(std::move(wp));
            }
        }

        // --- Step 5: Inverse transforms (U2†, U1) — quad phases, must be CPU ---
        for (auto& wp : tl_newWPs) {
            ApplyU2Dagger(wp, this->initTime, this->finalTime);
            ApplyU1(wp, this->finalTime, this->initTime);
            wp->SetTime(this->finalTime);
        }

        // Transfer tl_newWPs to atom ownership via buffer swap.
        // SwapWavePackets: destroys old WPs in atom, swaps buffers so the atom
        // holds the new WPs and tl_newWPs gets the old empty buffer back —
        // no heap allocation, and tl_newWPs retains its capacity for the next atom.
        atom->SwapWavePackets(tl_newWPs);
    }
}
