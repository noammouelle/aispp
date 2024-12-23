#ifndef AISWAVERONTS_HH
#define AISWAVERONTS_HH

#include <array>
#include <cmath>
#include "AISConstants.hh"

// flat wavefront and gradient wavefront
double flatWavefront(const std::array<double, 3>& pos, const double& param1, const double& param2);
std::array<double, 3> flatGradientWavefront(const std::array<double, 3>& pos, const double& param1, const double& param2);

// gaussian wavefront and gradient
double gaussianWavefront(const std::array<double, 3>& pos, const double& kz, const double& w0);
std::array<double, 3> gaussianGradientWavefront(const std::array<double, 3>& pos, const double& kz, const double& w0);

#endif