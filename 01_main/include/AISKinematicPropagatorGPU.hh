#ifndef AISKINEMATICPROPAGATOR_GPU_HH
#define AISKINEMATICPROPAGATOR_GPU_HH

#ifdef USE_CUDA

#include <vector>
#include "AISUtilities.hh"   // doubleThreeVector

// Launch the kinematic propagation CUDA kernel.
//
// Applicable only when the potential produces a CONSTANT acceleration
// (zero_pot or linear_pot / uniform gravity).  For position-dependent
// potentials (quadratic_pot) the caller must fall back to the CPU path.
//
// The kernel implements the exact analytic trajectory update:
//   pos_new = pos + vel*dt + 0.5*a*dt^2
//   vel_new = vel + a*dt
//
// Arrays are host pointers of length N (all wavepackets, all atoms, flat).
// pos/vel arrays are updated IN PLACE.
void launchKinematicKernel(
    double* pos_x, double* pos_y, double* pos_z,
    double* vel_x, double* vel_y, double* vel_z,
    double ax, double ay, double az,   // constant acceleration = -dUdx/m
    double dt,
    int N);

#endif // USE_CUDA
#endif // AISKINEMATICPROPAGATOR_GPU_HH
