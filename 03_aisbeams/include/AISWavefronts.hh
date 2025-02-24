#ifndef AISWAVERONTS_HH
#define AISWAVERONTS_HH

#include <array>
#include <cmath>
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

#endif