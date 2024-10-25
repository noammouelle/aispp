#ifndef AISPOTENTIAL_HH
#define AISPOTENTIAL_HH

#include <array>

// no gravity
double zeroU(std::array<double,3> pos, std::array<double,3> vel)
{
    return 0.0;
}
std::array<double,3> zeroGrad(std::array<double,3> pos, std::array<double,3> vel)
{
    return {0.,0.,0.};
}
std::array<std::array<double,3>,3> zeroHess(std::array<double,3> pos, std::array<double,3> vel)
{
    return {{{0.,0.,0.},
             {0.,0.,0.},
             {0.,0.,0.}}};
}

#endif