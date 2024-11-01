#include "AISWavefronts.hh"

// flat wavefront and gradient wavefront
double flatWavefront(const std::array<double, 3>& pos)
{
    return 0.0;
}
std::array<double, 3> flatGradientWavefront(const std::array<double, 3>& pos)
{
    return {0.0, 0.0, 0.0};
}