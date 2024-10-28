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