#include "AISPotentials.hh"

// no gravity
double zeroU(const std::array<double,3>& pos, const std::array<double,3>& vel)
{
    return 0.0;
}
std::array<double,3> zeroGrad(const std::array<double,3>& pos, const std::array<double,3>& vel)
{
    return {0.,0.,0.};
}
std::array<std::array<double,3>,3> zeroHess(const std::array<double,3>& pos, const std::array<double,3>& vel)
{
    return {{{0.,0.,0.},
             {0.,0.,0.},
             {0.,0.,0.}}};
}

// uniform gravity
double uniformGravityU(const std::array<double,3>& pos, const std::array<double,3>& vel)
{
    return g * pos[2];
}
std::array<double,3> uniformGravityGrad(const std::array<double,3>& pos, const std::array<double,3>& vel)
{
    return {0.,0.,massSr87 * g};
}

// linear gravity
double linearGravityU(const std::array<double,3>& pos, const std::array<double,3>& vel)
{
    return g * pos[2] + 0.5 * (g/radiusEarth * pos[0] * pos[0] + g/radiusEarth * pos[1] * pos[1] - 2 * g/radiusEarth * pos[2] * pos[2]);
}
std::array<double,3> linearGravityGrad(const std::array<double,3>& pos, const std::array<double,3>& vel)
{
    return {massSr87 * g/radiusEarth * pos[0], 
            massSr87 * g/radiusEarth * pos[1], 
            massSr87 * g - massSr87 * 2 * g/radiusEarth * pos[2]};
}
std::array<std::array<double,3>,3> linearGravityHess(const std::array<double,3>& pos, const std::array<double,3>& vel)
{
    return {{{massSr87 * g/radiusEarth, 0., 0.},
             {0., massSr87 * g/radiusEarth, 0.},
             {0., 0., -2 * massSr87 * g/radiusEarth}}};
}