#ifndef AISPOTENTIAL_HH
#define AISPOTENTIAL_HH

#include <array>

// no gravity
double zeroU(const std::array<double,3>& pos, const std::array<double,3>& vel);
std::array<double,3> zeroGrad(const std::array<double,3>& pos, const std::array<double,3>& vel);
std::array<std::array<double,3>,3> zeroHess(const std::array<double,3>& pos, const std::array<double,3>& vel);

#endif