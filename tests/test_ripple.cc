// wtype "ripple": uniform Rabi frequency with a prescribed sinusoidal transverse phase.
// Checks GetPhi, GetDelPhi against a central finite difference, GetRabiFreq, and the eps = 0 limit.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "AISLaserBeam.hh"
#include "AISUtilities.hh"

static int fails = 0;
static void check(bool ok, const char* what, double a, double b)
{
    if(!ok) { std::printf("FAIL %s: %.12g vs %.12g\n", what, a, b); ++fails; }
    else    { std::printf("ok   %s\n", what); }
}

int main()
{
    const double k = 8995953.76, q = 2.0*M_PI/4.0e-6, eps = 0.1, psi = 0.3, phi0 = 0.7, rabi = 2.0*M_PI*1000.0;
    AISLaserBeam beam({0.0, 0.0, k}, (__float128)2.696919e15, rabi, phi0);
    beam.SetW0(1.0e-3); beam.SetFocalLength(0.0); beam.SetZLaser(0.0); beam.SetBeamRadius(1.0);
    beam.SetBeamType("ripple");
    beam.SetRipple(q, eps, psi);
    const double xs[] = {-1.7e-6, 0.0, 0.53e-6, 2.9e-6};
    for(double x : xs)
    {
        doubleThreeVector p = {x, 0.2e-6, 3.0e-6};
        double phi = beam.GetPhi(p);
        check(std::fabs(phi - (phi0 + eps*std::cos(q*x + psi))) < 1e-12, "GetPhi analytic", phi, phi0 + eps*std::cos(q*x + psi));
        doubleThreeVector g = beam.GetDelPhi(p);
        const double h = 1.0e-9;
        doubleThreeVector pp = {x + h, p[1], p[2]}, pm = {x - h, p[1], p[2]};
        double fd = (beam.GetPhi(pp) - beam.GetPhi(pm))/(2.0*h);
        check(std::fabs(g[0] - fd) < 1e-6*std::fabs(eps*q) + 1e-9, "GetDelPhi x vs finite difference", g[0], fd);
        check(g[1] == 0.0 && g[2] == 0.0, "GetDelPhi y,z zero", g[1], g[2]);
        double om = beam.GetRabiFreq(p, (__float128)0, (__float128)1e-4);
        check(om == rabi, "GetRabiFreq uniform", om, rabi);
    }
    beam.SetRipple(q, 0.0, psi);
    doubleThreeVector p0 = {1.1e-6, 0.0, 0.0};
    check(beam.GetPhi(p0) == phi0, "eps = 0 gives phi0", beam.GetPhi(p0), phi0);
    check(beam.GetDelPhi(p0)[0] == 0.0, "eps = 0 gives zero gradient", beam.GetDelPhi(p0)[0], 0.0);
    std::printf("%s\n", fails ? "test_ripple: FAILED" : "test_ripple: all passed");
    return fails ? 1 : 0;
}
