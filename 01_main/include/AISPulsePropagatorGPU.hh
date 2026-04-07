#ifndef AISPULSEPROPAGATOR_GPU_HH
#define AISPULSEPROPAGATOR_GPU_HH

// This header is only compiled when USE_CUDA is defined.
// It declares the CUDA kernel launch helper called from
// AISPulsePropagator::PropagateEnsembleGPU().

#ifdef USE_CUDA

#include <vector>

// Plain-old-data struct passed to the CUDA kernel.
// All fields are precomputed on the CPU before the kernel launch so that
// the kernel contains zero __float128 arithmetic.
struct GPUBeamParams
{
    // k vector at time t1 (k0 + 0.5*kChirp*t1) — same for all atoms
    double kT1_x, kT1_y, kT1_z;
    // wavevector chirp
    double kChirp_x, kChirp_y, kChirp_z;
    // central Rabi frequency (rad/s)
    double rabiFreq;
    // Gaussian beam waist (m)
    double w0;
    // focal length (m)  — sign of k[2] applied separately
    double focalLength;
    // dOmega = double(omega - omegaSr87)  [computed via quad on CPU, once]
    double dOmegaScalar;
    // -omegaChirp * t1
    double omegaChirpTerm;
    // |kT1|^2 * hbar / (2*massSr87)  — recoil, constant since gradPhi==0
    double recoilTerm;
    // pulse window
    double t0, t1;
    // +1 or -1 depending on k[2] sign (for focal-length shift in Gaussian beam)
    int kzSign;
    // 1 = flat_square (Omega = rabiFreq everywhere)
    // 0 = gaussian
    int beamTypeFlat;
};

// Launch the U3-UltraFast CUDA kernel.
// Inputs  (length N): per-wavepacket position, velocity, state
// Outputs (length N): amplitude and phase-delta for both daughter wavepackets
void launchU3UltraFastKernel(
    const double* pos_x, const double* pos_y, const double* pos_z,
    const double* vel_x, const double* vel_y, const double* vel_z,
    const double* amp_in,
    const int*    state_in,
    double* amp0_out,    double* dphase0_out,   // "keep-state" wavepacket
    double* amp1_out,    double* dphase1_out,   // "flip-state" wavepacket
    const GPUBeamParams& bp,
    int N);

#endif // USE_CUDA
#endif // AISPULSEPROPAGATOR_GPU_HH
