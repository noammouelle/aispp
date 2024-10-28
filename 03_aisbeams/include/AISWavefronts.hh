#ifndef AISWAVERONTS_HH
#define AISWAVERONTS_HH

#include <array>

// flat wavefront and gradient wavefront
double flatWavefront(std::array<double, 3> pos);
std::array<double, 3> flatGradientWavefront(std::array<double, 3> pos);

#endif