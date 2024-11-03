#ifndef AISPOTENTIAL_HH
#define AISPOTENTIAL_HH

#include <array>
#include <iostream>

#include "AISConstants.hh"

// no gravity
double zeroU(const std::array<double,3>& pos, const std::array<double,3>& vel);
std::array<double,3> zeroGrad(const std::array<double,3>& pos, const std::array<double,3>& vel);
std::array<std::array<double,3>,3> zeroHess(const std::array<double,3>& pos, const std::array<double,3>& vel);

// uniform gravity
double uniformGravityU(const std::array<double,3>& pos, const std::array<double,3>& vel);
std::array<double,3> uniformGravityGrad(const std::array<double,3>& pos, const std::array<double,3>& vel);

// linear gravity
double linearGravityU(const std::array<double,3>& pos, const std::array<double,3>& vel);
std::array<double,3> linearGravityGrad(const std::array<double,3>& pos, const std::array<double,3>& vel);
std::array<std::array<double,3>,3> linearGravityHess(const std::array<double,3>& pos, const std::array<double,3>& vel);

#endif