#ifndef AISPOTENTIAL_HH
#define AISPOTENTIAL_HH

#include <array>
#include <iostream>

#include "AISConstants.hh"
#include "AISUtilities.hh"

// Convention reminder for everything in this file:
//   * U(...)        returns the potential energy DIVIDED BY the atomic mass.
//   * dUdx(...)     returns the full gradient dU/dx (not divided by mass).
//   * dUdp(...)     returns dU/dp, which has units of velocity.
//   * d2Udxdx, d2Udxdp, d2Udpdp are the corresponding second derivatives,
//     always with respect to the canonical momentum p (never to v = p/m).
//   * The second argument of every function is the wavepacket's canonical
//     velocity v = p/m, which equals the frame velocity xdot only when the
//     potential has no velocity dependence (i.e. when dUdp == 0).

// no gravity
inline double zeroU(const std::array<double,3>& pos, const std::array<double,3>& vel)
{
    return 0.0;
}
inline std::array<double,3> zeroGrad(const std::array<double,3>& pos, const std::array<double,3>& vel)
{
    return {0.,0.,0.};
}
inline std::array<std::array<double,3>,3> zeroHess(const std::array<double,3>& pos, const std::array<double,3>& vel)
{
    return {{{0.,0.,0.},
             {0.,0.,0.},
             {0.,0.,0.}}};
}

// uniform gravity
inline double uniformGravityU(const std::array<double,3>& pos, const std::array<double,3>& vel)
{
    return g * pos[2];
}
inline std::array<double,3> uniformGravityGrad(const std::array<double,3>& pos, const std::array<double,3>& vel)
{
    return {0.,0.,massSr87 * g};
}

// linear gravity
inline double linearGravityU(const std::array<double,3>& pos, const std::array<double,3>& vel)
{
    return g * pos[2] + 0.5 * (g/radiusEarth * pos[0] * pos[0] + g/radiusEarth * pos[1] * pos[1] - 2 * g/radiusEarth * pos[2] * pos[2]);
}
inline std::array<double,3> linearGravityGrad(const std::array<double,3>& pos, const std::array<double,3>& vel)
{
    return {massSr87 * g/radiusEarth * pos[0],
            massSr87 * g/radiusEarth * pos[1],
            massSr87 * g - massSr87 * 2 * g/radiusEarth * pos[2]};
}
inline std::array<std::array<double,3>,3> linearGravityHess(const std::array<double,3>& pos, const std::array<double,3>& vel)
{
    return {{{massSr87 * g/radiusEarth, 0., 0.},
             {0., massSr87 * g/radiusEarth, 0.},
             {0., 0., -2 * massSr87 * g/radiusEarth}}};
}

// ─── Rotating frame (Coriolis + centrifugal) ──────────────────────────────────
//
// In a frame rotating at constant angular velocity Omega the Lagrangian is
//
//     L = (1/2) m |xdot + Omega x x|^2 - V(x),
//
// so the canonical momentum is p = m (xdot + Omega x x) and the Hamiltonian is
//
//     H = p^2/(2m) - Omega.(x x p) + V(x).
//
// AIS++ integrates Hamilton's equations in the form
//     xdot = v + dU/dp,      vdot = -(dU/dx)/m,      with v = p/m,
// so the rotating frame enters entirely through the velocity-dependent piece of
// U = H - p^2/(2m):
//
//     U(x,v)/m = V(x)/m - Omega.(x x v)
//     dU/dx    = grad V + m (Omega x v)
//     dU/dp    = -(Omega x x)
//     d2U/dxdp = [Omega]x            (skew-symmetric cross-product matrix)
//     d2U/dpdp = 0
//
// The centrifugal force is NOT added here. It is generated automatically by
// p^2/(2m) once p is the canonical momentum; adding it explicitly would double
// count it. Substituting the above into xdot/vdot and eliminating p recovers
//     m xddot = -grad V - 2 m (Omega x xdot) - m Omega x (Omega x x).
//
// Because U is now velocity dependent, v = p/m is NOT the frame velocity. The
// frame velocity is xdot = v - Omega x x. Conversions live at the I/O boundary
// (see frameFromCanonicalVelocity / canonicalFromFrameVelocity below); the
// physics core works exclusively in canonical variables.

// Angular velocity of the rotating frame, in rad/s, expressed in the simulation
// frame. Defined in AISPotentials.cc. The potential functions are plain function
// pointers and cannot capture state, so this is set once at start-up by
// AISDriver via SetRotationRate().
extern doubleThreeVector gRotationRate;

void SetRotationRate(const doubleThreeVector& omega);

// Frame velocity xdot = v - Omega x x, given the canonical velocity v = p/m.
inline doubleThreeVector frameFromCanonicalVelocity(const doubleThreeVector& pos,
                                                    const doubleThreeVector& velCanonical)
{
    return matrixAdd(velCanonical, scalarMultiply(crossProduct(gRotationRate, pos), -1.0));
}

// Canonical velocity v = p/m = xdot + Omega x x, given the frame velocity xdot.
inline doubleThreeVector canonicalFromFrameVelocity(const doubleThreeVector& pos,
                                                    const doubleThreeVector& velFrame)
{
    return matrixAdd(velFrame, crossProduct(gRotationRate, pos));
}

// -- rotation-only pieces, shared by every rotating_* potential ----------------

// -Omega.(x x v)   (already per unit mass)
inline double rotationU(const std::array<double,3>& pos, const std::array<double,3>& vel)
{
    return -dotProduct(gRotationRate, crossProduct(pos, vel));
}
// +m (Omega x v)
inline std::array<double,3> rotationGrad(const std::array<double,3>& pos, const std::array<double,3>& vel)
{
    return scalarMultiply(crossProduct(gRotationRate, vel), massSr87);
}
// -(Omega x x)
inline std::array<double,3> rotationDUdp(const std::array<double,3>& pos, const std::array<double,3>& vel)
{
    return scalarMultiply(crossProduct(gRotationRate, pos), -1.0);
}
// [Omega]x
inline std::array<std::array<double,3>,3> rotationD2Udxdp(const std::array<double,3>& pos, const std::array<double,3>& vel)
{
    return crossMatrix(gRotationRate);
}

// -- rotation only (no gravity) -----------------------------------------------
inline double rotatingU(const std::array<double,3>& pos, const std::array<double,3>& vel)
{
    return rotationU(pos, vel);
}
inline std::array<double,3> rotatingGrad(const std::array<double,3>& pos, const std::array<double,3>& vel)
{
    return rotationGrad(pos, vel);
}

// -- rotation + uniform gravity -----------------------------------------------
inline double rotatingUniformGravityU(const std::array<double,3>& pos, const std::array<double,3>& vel)
{
    return uniformGravityU(pos, vel) + rotationU(pos, vel);
}
inline std::array<double,3> rotatingUniformGravityGrad(const std::array<double,3>& pos, const std::array<double,3>& vel)
{
    return matrixAdd(uniformGravityGrad(pos, vel), rotationGrad(pos, vel));
}

// -- rotation + gravity gradient ----------------------------------------------
inline double rotatingLinearGravityU(const std::array<double,3>& pos, const std::array<double,3>& vel)
{
    return linearGravityU(pos, vel) + rotationU(pos, vel);
}
inline std::array<double,3> rotatingLinearGravityGrad(const std::array<double,3>& pos, const std::array<double,3>& vel)
{
    return matrixAdd(linearGravityGrad(pos, vel), rotationGrad(pos, vel));
}

#endif
