#include "AISWavefronts.hh"

// flat wavefront and gradient wavefront
double flatWavefront(std::array<double, 3> pos)
{
    return 0.0;
}
std::array<double, 3> flatGradientWavefront(std::array<double, 3> pos)
{
    return {0.0, 0.0, 0.0};
}