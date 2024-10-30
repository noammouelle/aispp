#ifndef AISWAVERONTS_HH
#define AISWAVERONTS_HH

#include <array>

// flat wavefront and gradient wavefront
double flatWavefront(const std::array<double, 3>& pos);
std::array<double, 3> flatGradientWavefront(const std::array<double, 3>& pos);

#endif