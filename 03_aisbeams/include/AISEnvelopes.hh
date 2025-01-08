#ifndef AISENVELOPES_HH
#define AISENVELOPES_HH

#include <array>
#include <cmath>
#include <iostream>

#include "AISConstants.hh"

// square pulse and flat envelope
double flatSquareEnvelope(const std::array<double, 3>& pos, 
                          const double& param1,
                          const __float128& t0, const __float128& t1);

// gaussian envelope
double gaussianEnvelope(const std::array<double, 3>& pos, 
                        const double& w0,
                        const __float128& t0, const __float128& t1);

#endif