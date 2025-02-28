#ifndef AISPOTENTIAL_HH
#define AISPOTENTIAL_HH

#include <array>
#include <iostream>

#include "AISConstants.hh"

// no gravity
inline double zeroU(const std::array<double,3>& pos, const std::array<double,3>& vel)
{
    return 0.0;
}
inline std::array<double,3> zeroGrad(const std::array<double,3>& pos, const std::array<double,3>& vel)
{
    return {0.,0.,0.};
}
inline std::array<std::array<double,3>,3> zeroHess(const std::array<double,3>& pos, const std::array<double,3>& vel)
{
    return {{{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}}};
}

// uniform gravity
inline double uniformGravityU(const std::array<double,3>& pos, const std::array<double,3>& vel)
{
    return g * pos[2];
}
inline std::array<double,3> uniformGravityGrad(const std::array<double,3>& pos, const std::array<double,3>& vel)
{
    return {0.,0.,massSr87 * g};
}

// linear gravity
double linearGravityU(const std::array<double,3>& pos, const std::array<double,3>& vel);
std::array<double,3> linearGravityGrad(const std::array<double,3>& pos, const std::array<double,3>& vel);
inline std::array<std::array<double,3>,3> linearGravityHess(const std::array<double,3>& pos, const std::array<double,3>& vel)
{
    return {{{massSr87 * g/radiusEarth, 0., 0.},
             {0., massSr87 * g/radiusEarth, 0.},
             {0., 0., -2 * massSr87 * g/radiusEarth}}};
}

#endif