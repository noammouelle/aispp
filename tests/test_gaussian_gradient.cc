// GetDelPhi for wtype "gaussian" must equal the finite-difference gradient of GetPhi (both beam directions, at and off the waist).
#include <cmath>
#include <cstdio>
#include "AISLaserBeam.hh"
#include "AISUtilities.hh"
int main()
{
    int fails = 0; const double k = 8995953.76, w0 = 20e-6, zR = M_PI*w0*w0/698.4457e-9;
    for (double kz : {k, -k}) for (double f : {0.0, zR, -2.5*zR}) {
        AISLaserBeam b({0.0, 0.0, kz}, (__float128)2.69e15, 2*M_PI*1000.0, 0.3);
        b.SetW0(w0); b.SetFocalLength(f); b.SetZLaser(0.0); b.SetBeamRadius(1.0); b.SetBeamType("gaussian");
        for (double x : {0.0, 7e-6, -19e-6}) for (double z : {0.0, 3e-6}) {
            doubleThreeVector p = {x, 4e-6, z}; doubleThreeVector g = b.GetDelPhi(p); double worst = 0;
            for (int a = 0; a < 3; ++a) { const double h = 1e-9; doubleThreeVector pp = p, pm = p; pp[a] += h; pm[a] -= h;
                double fd = (b.GetPhi(pp) - b.GetPhi(pm))/(2*h); worst = fmax(worst, fabs(g[a] - fd)/(fabs(fd) + 1.0)); }
            if (worst > 1e-5) { ++fails; std::printf("FAIL kz %+.0f f %.3g x %.1e z %.1e: rel err %.2e\n", kz, f, x, z, worst); }
        }
    }
    std::printf("%s\n", fails ? "test_gaussian_gradient: FAILED" : "test_gaussian_gradient: all passed"); return fails ? 1 : 0;
}
