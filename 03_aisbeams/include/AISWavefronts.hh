#ifndef AISWAVERONTS_HH
#define AISWAVERONTS_HH

#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include "AISConstants.hh"
#include "AISUtilities.hh"

// flat wavefront and gradient wavefront
double flatWavefront(const std::array<double, 3>& pos, const double& param1, const double& param2);
std::array<double, 3> flatGradientWavefront(const std::array<double, 3>& pos, const double& param1, const double& param2);

// gaussian wavefront and gradient
double gaussianWavefront(const std::array<double, 3>& pos, const double& kz, const double& w0);
std::array<double, 3> gaussianGradientWavefront(const std::array<double, 3>& pos, const double& kz, const double& w0);

// utility functions to compute the zernike polynomial type aberrations
double Rmn(int m, int n, double rho);
double Zmn(int m, int n, double rho, double theta);
double a(int n, int m, int nprime, int mprime);
std::array<double, 3> GradZmn(int m, int n, double rho, double theta); // cartesian gradient of the zernike polynomial

double sgn(double val);

// ---------------------------------------------------------------------------
// Analytic single transverse Fourier mode on a Gaussian carrier.
//
// A mirror surface ripple imprints phi = eps cos(q.rho + psi) at the mirror
// plane. For small eps the reflected field is a carrier plus two sidebands,
// each of which is the carrier times a plane wave -- a tilted Gaussian, which
// propagates in closed form. Over a distance dz from the mirror:
//
//   E(rho,z) = E0(rho,z)
//            + (i eps/2) exp(-i q^2 dz / 2k)
//              [ exp(+i(q.rho+psi)) E0(rho - q dz/k, z)
//              + exp(-i(q.rho+psi)) E0(rho + q dz/k, z) ]
//
// Each sideband walks off transversely by q dz / k and picks up q^2 dz / 2k,
// a phase that advances by 2*pi over one Talbot length 2 lambda_perp^2/lambda.
// The Talbot effect is therefore exact here rather than something a sampled
// z grid has to resolve, and there is no grid to leave: the field is defined
// at every point, for a source at 0-20 m or at 100-120 m alike.
//
// Valid to first order in eps; the neglected term is smaller by ~eps/4
// (4.5e-3 for a 1 nm ripple at 698 nm).
// ---------------------------------------------------------------------------

struct FourierMode
{
    double qx      = 0.0;  // rad/m
    double qy      = 0.0;  // rad/m
    double epsilon = 0.0;  // peak optical phase amplitude at the mirror, rad
    double psi     = 0.0;  // spatial phase of the ripple, rad
    double mirrorZ = 0.0;  // plane the ripple is imprinted at, m

    bool active() const { return epsilon != 0.0 && (qx != 0.0 || qy != 0.0); }
    double q() const { return std::sqrt(qx*qx + qy*qy); }
};

// Complex carrier: gaussianEnvelope * exp(i * gaussianWavefront), in the same
// shifted-coordinate convention the caller uses for the unaberrated beam.
std::complex<double> gaussianCarrierField(const std::array<double, 3>& shiftedPos,
                                          const double& kz, const double& w0);

// Complex field of carrier + mode. shiftedPos is the carrier coordinate;
// transversePos carries the LAB transverse coordinates the ripple is a
// function of, and dzFromMirror the signed distance from the mirror plane.
std::complex<double> fourierModeField(const std::array<double, 3>& shiftedPos,
                                      const double& transverseX,
                                      const double& transverseY,
                                      const double& dzFromMirror,
                                      const double& kz, const double& w0,
                                      const FourierMode& mode);

#endif