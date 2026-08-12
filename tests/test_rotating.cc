// Unit tests for the rotating-frame (Coriolis) potentials and for the
// velocity-dependent terms they switch on in AISKinematicPropagator.
//
// Build:  cmake --build build --target test_rotating
// Run:    ./build/tests/test_rotating
//
// Everything here is analytic or finite-difference: no input files, no HDF5.

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "AISPotentials.hh"
#include "AISKinematicPropagator.hh"

// ── tiny test harness ────────────────────────────────────────────────────────

static int gChecks = 0;
static int gFailures = 0;

static void check(bool ok, const std::string& what)
{
    ++gChecks;
    if (!ok) { ++gFailures; std::printf("  FAIL  %s\n", what.c_str()); }
}

static void checkClose(double got, double want, double tol, const std::string& what)
{
    ++gChecks;
    double err = std::fabs(got - want);
    double scale = std::fmax(1.0, std::fmax(std::fabs(got), std::fabs(want)));
    if (err > tol * scale)
    {
        ++gFailures;
        std::printf("  FAIL  %s: got %.17g want %.17g (rel err %.3g > %.3g)\n",
                    what.c_str(), got, want, err / scale, tol);
    }
}

static void checkVecClose(const doubleThreeVector& got, const doubleThreeVector& want,
                          double tol, const std::string& what)
{
    for (int i = 0; i < 3; ++i)
        checkClose(got[i], want[i], tol, what + "[" + std::to_string(i) + "]");
}

static void section(const std::string& name)
{
    std::printf("\n%s\n", name.c_str());
}

// ── helpers ──────────────────────────────────────────────────────────────────

// Rodrigues rotation of v about unit axis n by angle theta.
static doubleThreeVector rotateAbout(const doubleThreeVector& v,
                                     const doubleThreeVector& omega, double t)
{
    double norm = std::sqrt(dotProduct(omega, omega));
    if (norm == 0.0) return v;
    doubleThreeVector n = scalarMultiply(omega, 1.0 / norm);
    double theta = norm * t;
    double c = std::cos(theta), s = std::sin(theta);
    doubleThreeVector term1 = scalarMultiply(v, c);
    doubleThreeVector term2 = scalarMultiply(crossProduct(n, v), s);
    doubleThreeVector term3 = scalarMultiply(n, dotProduct(n, v) * (1.0 - c));
    return matrixAdd(term1, matrixAdd(term2, term3));
}

using PotFn  = double(*)(const doubleThreeVector&, const doubleThreeVector&);
using GradFn = doubleThreeVector(*)(const doubleThreeVector&, const doubleThreeVector&);
using HessFn = double3x3Matrix(*)(const doubleThreeVector&, const doubleThreeVector&);

static AISKinematicPropagator makePropagator(PotFn U, GradFn dUdx, GradFn dUdp,
                                             HessFn d2Udxdx, HessFn d2Udxdp, HessFn d2Udpdp)
{
    AISKinematicPropagator prop(std::make_shared<PotFn>(U),
                                std::make_shared<GradFn>(dUdx),
                                std::make_shared<GradFn>(dUdp),
                                std::make_shared<HessFn>(d2Udxdx),
                                std::make_shared<HessFn>(d2Udxdp),
                                std::make_shared<HessFn>(d2Udpdp));
    prop.ultraFast = false;
    prop.odeAbsTol = 1e-13;
    prop.odeRelTol = 1e-13;
    prop.qagAbsTol = 1e-13;
    prop.qagRelTol = 1e-13;
    return prop;
}

// ── 1. derivative consistency by finite differences ──────────────────────────
//
// Conventions (see AISPotentials.hh): U returns U/m, dUdx is the full gradient,
// dUdp = dU/dp has units of velocity, and the second derivatives are taken with
// respect to p, so d/dp = (1/m) d/dv.

static void testDerivatives(PotFn U, GradFn dUdx, GradFn dUdp,
                            HessFn d2Udxdx, HessFn d2Udxdp, HessFn d2Udpdp,
                            const std::string& name)
{
    section("Derivative consistency (finite differences): " + name);

    doubleThreeVector x = {0.031, -0.017, 0.24};
    doubleThreeVector v = {0.013, 0.0071, -3.1};
    const double hx = 1e-6, hv = 1e-6;

    for (int j = 0; j < 3; ++j)
    {
        doubleThreeVector xp = x, xm = x;
        xp[j] += hx; xm[j] -= hx;
        doubleThreeVector vp = v, vm = v;
        vp[j] += hv; vm[j] -= hv;

        // dU/dx_j = m * d(U/m)/dx_j
        checkClose(dUdx(x, v)[j], massSr87 * (U(xp, v) - U(xm, v)) / (2 * hx), 1e-6,
                   name + " dUdx[" + std::to_string(j) + "]");

        // dU/dp_j = d(U/m)/dv_j
        checkClose(dUdp(x, v)[j], (U(x, vp) - U(x, vm)) / (2 * hv), 1e-6,
                   name + " dUdp[" + std::to_string(j) + "]");

        for (int i = 0; i < 3; ++i)
        {
            checkClose(d2Udxdx(x, v)[i][j], (dUdx(xp, v)[i] - dUdx(xm, v)[i]) / (2 * hx), 1e-5,
                       name + " d2Udxdx[" + std::to_string(i) + "][" + std::to_string(j) + "]");

            // d2U/dx_i dp_j = (1/m) d(dUdx_i)/dv_j
            checkClose(d2Udxdp(x, v)[i][j],
                       (dUdx(x, vp)[i] - dUdx(x, vm)[i]) / (2 * hv) / massSr87, 1e-5,
                       name + " d2Udxdp[" + std::to_string(i) + "][" + std::to_string(j) + "]");

            // d2U/dp_i dp_j = (1/m) d(dUdp_i)/dv_j
            checkClose(d2Udpdp(x, v)[i][j],
                       (dUdp(x, vp)[i] - dUdp(x, vm)[i]) / (2 * hv) / massSr87, 1e-5,
                       name + " d2Udpdp[" + std::to_string(i) + "][" + std::to_string(j) + "]");

            // mixed partials must commute: d2U/dp dx = transpose(d2U/dx dp)
            checkClose((dUdp(xp, v)[i] - dUdp(xm, v)[i]) / (2 * hx),
                       d2Udxdp(x, v)[j][i], 1e-5,
                       name + " d2Udpdx symmetry [" + std::to_string(i) + "][" + std::to_string(j) + "]");
        }
    }
}

// ── 2. free particle in a rotating frame vs the analytic solution ────────────
//
// With V = 0 the inertial motion is a straight line. The canonical velocity
// v = p/m is the inertial velocity written in rotating-frame components, so
//     r_rot(t) = R(-Omega t) [r0 + v_canon(0) t],
//     v_canon(t) = R(-Omega t) v_canon(0),
// where R(theta) rotates about Omega. This also confirms that the frame
// velocity xdot = v - Omega x r is the thing that starts at the requested v0.

static void testFreeParticleTrajectory()
{
    section("Free particle in a rotating frame vs analytic solution");

    doubleThreeVector omega = {0.0, 7.292115e-5, 0.0};   // Earth rate about y
    SetRotationRate(omega);

    auto prop = makePropagator(rotatingU, rotatingGrad, rotationDUdp,
                               zeroHess, rotationD2Udxdp, zeroHess);

    doubleThreeVector x0 = {0.004, -0.002, 0.11};
    doubleThreeVector v0Frame = {0.013, -0.006, 2.7};
    doubleThreeVector v0Canon = canonicalFromFrameVelocity(x0, v0Frame);

    // the conversion must round-trip
    checkVecClose(frameFromCanonicalVelocity(x0, v0Canon), v0Frame, 1e-12,
                  "frame/canonical velocity round trip");

    const double t1 = 2.2;
    auto res = prop.CalculateNewPhaseSpaceCoords(0.0, t1, x0, v0Canon);

    doubleThreeVector rInertial = matrixAdd(x0, scalarMultiply(v0Canon, t1));
    doubleThreeVector wantPos = rotateAbout(rInertial, omega, -t1);
    doubleThreeVector wantVel = rotateAbout(v0Canon, omega, -t1);

    checkVecClose(res[0], wantPos, 1e-9, "rotating free-particle position");
    checkVecClose(res[1], wantVel, 1e-9, "rotating free-particle canonical velocity");

    // |v_canon| is conserved (it is the inertial speed)
    checkClose(std::sqrt(dotProduct(res[1], res[1])),
               std::sqrt(dotProduct(v0Canon, v0Canon)), 1e-11,
               "canonical speed conserved");

    SetRotationRate({0.0, 0.0, 0.0});
}

// ── 3. equations of motion vs an independent RK4 in physical variables ───────
//
// Integrates m xddot = -grad V - 2 m (Omega x xdot) - m Omega x (Omega x x)
// directly, in frame variables, with a completely separate integrator. This is
// an independent implementation of the physics, not a restatement of the code.

static void testAgainstIndependentIntegrator()
{
    section("Rotating frame + gravity vs independent RK4 in frame variables");

    doubleThreeVector omega = {1.3e-4, -7.0e-5, 4.1e-5};
    SetRotationRate(omega);

    auto prop = makePropagator(rotatingUniformGravityU, rotatingUniformGravityGrad,
                               rotationDUdp, zeroHess, rotationD2Udxdp, zeroHess);

    doubleThreeVector x0 = {0.003, 0.001, 0.0};
    doubleThreeVector v0Frame = {0.01, -0.02, 3.0};

    // reference: RK4 on (x, xdot) with the textbook rotating-frame acceleration
    auto accel = [&](const doubleThreeVector& x, const doubleThreeVector& xd) {
        doubleThreeVector coriolis = scalarMultiply(crossProduct(omega, xd), -2.0);
        doubleThreeVector centrifugal =
            scalarMultiply(crossProduct(omega, crossProduct(omega, x)), -1.0);
        doubleThreeVector gravity = {0.0, 0.0, -g};
        return matrixAdd(gravity, matrixAdd(coriolis, centrifugal));
    };

    const double T = 2.0;
    const int nSteps = 400000;
    const double dt = T / nSteps;
    doubleThreeVector x = x0, xd = v0Frame;
    for (int s = 0; s < nSteps; ++s)
    {
        doubleThreeVector k1x = xd,                                  k1v = accel(x, xd);
        doubleThreeVector k2x = matrixAdd(xd, scalarMultiply(k1v, dt/2)),
                          k2v = accel(matrixAdd(x, scalarMultiply(k1x, dt/2)),
                                      matrixAdd(xd, scalarMultiply(k1v, dt/2)));
        doubleThreeVector k3x = matrixAdd(xd, scalarMultiply(k2v, dt/2)),
                          k3v = accel(matrixAdd(x, scalarMultiply(k2x, dt/2)),
                                      matrixAdd(xd, scalarMultiply(k2v, dt/2)));
        doubleThreeVector k4x = matrixAdd(xd, scalarMultiply(k3v, dt)),
                          k4v = accel(matrixAdd(x, scalarMultiply(k3x, dt)),
                                      matrixAdd(xd, scalarMultiply(k3v, dt)));
        doubleThreeVector dx = scalarMultiply(matrixAdd(matrixAdd(k1x, scalarMultiply(k2x, 2.0)),
                                                        matrixAdd(scalarMultiply(k3x, 2.0), k4x)), dt/6);
        doubleThreeVector dv = scalarMultiply(matrixAdd(matrixAdd(k1v, scalarMultiply(k2v, 2.0)),
                                                        matrixAdd(scalarMultiply(k3v, 2.0), k4v)), dt/6);
        x = matrixAdd(x, dx);
        xd = matrixAdd(xd, dv);
    }

    auto res = prop.CalculateNewPhaseSpaceCoords(0.0, T, x0, canonicalFromFrameVelocity(x0, v0Frame));
    doubleThreeVector gotPos = res[0];
    doubleThreeVector gotFrameVel = frameFromCanonicalVelocity(gotPos, res[1]);

    // absolute tolerances: displacements are ~metres, the rotation correction ~mm
    for (int i = 0; i < 3; ++i)
    {
        checkClose(gotPos[i] - x[i], 0.0, 1e-9, "position residual vs RK4 [" + std::to_string(i) + "]");
        checkClose(gotFrameVel[i] - xd[i], 0.0, 1e-9, "frame velocity residual vs RK4 [" + std::to_string(i) + "]");
    }

    SetRotationRate({0.0, 0.0, 0.0});
}

// ── 4. the action in a rotating frame ────────────────────────────────────────
//
// In canonical variables the rotating-frame Lagrangian is exactly
//     L/m = v^2/2 - V/m,
// because the -Omega.(x x v) in U cancels against the v.(dU/dp) term of
// L = p.xdot - H. With V = 0 the canonical speed is conserved, so
//     dScl = (|v|^2 - |vTilde|^2) (t1 - t0) / 2
// in closed form. This is the test that fails if get_dL drops the velocity-
// dependent terms, and it is the term that carries the Sagnac phase.

static void testRotatingAction()
{
    section("Action integral in a rotating frame");

    doubleThreeVector omega = {0.0, 7.292115e-5, 0.0};
    SetRotationRate(omega);

    auto prop = makePropagator(rotatingU, rotatingGrad, rotationDUdp,
                               zeroHess, rotationD2Udxdp, zeroHess);
    prop.SetAddEnergyPhase(false);

    // the reference ("tilde") trajectory the propagator expands about
    doubleThreeVector atomPos0 = {0.0, 0.0, 0.0};
    doubleThreeVector atomVel0 = canonicalFromFrameVelocity(atomPos0, {0.0, 0.0, 3.0});

    // the wavepacket, displaced in both position and velocity
    doubleThreeVector pos0 = {0.002, -0.001, 0.0};
    doubleThreeVector vel0 = canonicalFromFrameVelocity(pos0, {0.011, 0.004, 3.05});

    const double t0 = 0.0, t1 = 1.5;
    std::array<double,2> scl = prop.get_dScl(t0, t1, pos0, vel0, atomPos0, atomVel0, 0.0q);

    double want = 0.5 * (dotProduct(vel0, vel0) - dotProduct(atomVel0, atomVel0)) * (t1 - t0);
    checkClose(scl[0], want, 1e-9, "rotating-frame action equals (v^2 - vTilde^2)(t1-t0)/2");

    SetRotationRate({0.0, 0.0, 0.0});
}

// ── 5. Omega = 0 must reproduce the inertial potentials bit for bit ──────────

static void testZeroRotationReduction()
{
    section("Omega = 0 reduces the rotating potentials to the inertial ones");

    SetRotationRate({0.0, 0.0, 0.0});

    doubleThreeVector x = {0.03, -0.02, 0.5};
    doubleThreeVector v = {0.1, 0.2, -3.0};

    check(rotatingU(x, v) == zeroU(x, v), "rotatingU == zeroU");
    check(rotatingUniformGravityU(x, v) == uniformGravityU(x, v),
          "rotatingUniformGravityU == uniformGravityU");
    check(rotatingLinearGravityU(x, v) == linearGravityU(x, v),
          "rotatingLinearGravityU == linearGravityU");

    for (int i = 0; i < 3; ++i)
    {
        check(rotatingGrad(x, v)[i] == zeroGrad(x, v)[i], "rotatingGrad == zeroGrad");
        check(rotatingUniformGravityGrad(x, v)[i] == uniformGravityGrad(x, v)[i],
              "rotatingUniformGravityGrad == uniformGravityGrad");
        check(rotatingLinearGravityGrad(x, v)[i] == linearGravityGrad(x, v)[i],
              "rotatingLinearGravityGrad == linearGravityGrad");
        check(rotationDUdp(x, v)[i] == 0.0, "rotationDUdp == 0");
        for (int j = 0; j < 3; ++j)
            check(rotationD2Udxdp(x, v)[i][j] == 0.0, "rotationD2Udxdp == 0");
    }

    // and the propagated trajectory must match the inertial one exactly
    auto rotProp = makePropagator(rotatingUniformGravityU, rotatingUniformGravityGrad,
                                  rotationDUdp, zeroHess, rotationD2Udxdp, zeroHess);
    auto inertialProp = makePropagator(uniformGravityU, uniformGravityGrad,
                                       zeroGrad, zeroHess, zeroHess, zeroHess);
    doubleThreeVector x0 = {0.001, 0.002, 0.0};
    doubleThreeVector v0 = {0.01, -0.01, 2.5};
    auto a = rotProp.CalculateNewPhaseSpaceCoords(0.0, 1.7, x0, v0);
    auto b = inertialProp.CalculateNewPhaseSpaceCoords(0.0, 1.7, x0, v0);
    checkVecClose(a[0], b[0], 1e-14, "Omega=0 position matches inertial");
    checkVecClose(a[1], b[1], 1e-14, "Omega=0 velocity matches inertial");

    // ... including the action
    auto sa = rotProp.get_dScl(0.0, 1.7, x0, v0, x0, v0, 0.0q);
    auto sb = inertialProp.get_dScl(0.0, 1.7, x0, v0, x0, v0, 0.0q);
    checkClose(sa[0], sb[0], 1e-14, "Omega=0 action matches inertial");
}

// ── 6. Coriolis deflection has the expected size and direction ───────────────
//
// A particle launched vertically for time T in a frame rotating about +y is
// deflected along x by -Omega_y * v_z * T^2 to leading order (the standard
// "eastward deflection" result, sign fixed by the frame orientation).

static void testCoriolisDeflectionScaling()
{
    section("Coriolis deflection magnitude and scaling");

    const double omegaY = 7.292115e-5;
    SetRotationRate({0.0, omegaY, 0.0});

    auto prop = makePropagator(rotatingU, rotatingGrad, rotationDUdp,
                               zeroHess, rotationD2Udxdp, zeroHess);

    doubleThreeVector x0 = {0.0, 0.0, 0.0};
    const double vz = 3.0;

    for (double T : {0.5, 1.0, 2.0})
    {
        doubleThreeVector v0 = canonicalFromFrameVelocity(x0, {0.0, 0.0, vz});
        auto res = prop.CalculateNewPhaseSpaceCoords(0.0, T, x0, v0);
        // leading order: x(T) = -Omega_y v_z T^2, exact to O((Omega T)^3)
        double want = -omegaY * vz * T * T;
        checkClose(res[0][0], want, 2e-4, "Coriolis deflection at T=" + std::to_string(T));
    }

    SetRotationRate({0.0, 0.0, 0.0});
}

// ── main ─────────────────────────────────────────────────────────────────────

int main()
{
    std::printf("ais++ rotating-frame unit tests\n");

    // finite-difference checks are done at a large, physically silly rotation
    // rate so the Coriolis terms are not lost in round-off
    SetRotationRate({0.7, -0.4, 1.1});
    testDerivatives(rotatingU, rotatingGrad, rotationDUdp,
                    zeroHess, rotationD2Udxdp, zeroHess, "rotating_pot");
    testDerivatives(rotatingUniformGravityU, rotatingUniformGravityGrad, rotationDUdp,
                    zeroHess, rotationD2Udxdp, zeroHess, "rotating_linear_pot");
    testDerivatives(rotatingLinearGravityU, rotatingLinearGravityGrad, rotationDUdp,
                    linearGravityHess, rotationD2Udxdp, zeroHess, "rotating_quadratic_pot");
    SetRotationRate({0.0, 0.0, 0.0});

    testFreeParticleTrajectory();
    testAgainstIndependentIntegrator();
    testRotatingAction();
    testZeroRotationReduction();
    testCoriolisDeflectionScaling();

    std::printf("\n%d checks, %d failures\n", gChecks, gFailures);
    if (gFailures == 0) std::printf("ALL TESTS PASSED\n");
    return gFailures == 0 ? 0 : 1;
}
