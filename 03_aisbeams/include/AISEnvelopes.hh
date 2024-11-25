#ifndef AISENVELOPES_HH
#define AISENVELOPES_HH

#include <array>
#include <cmath>

// square pulse and flat envelope
double flatSquareEnvelope(const std::array<double, 3>& pos, const __float128& t0, const __float128& t1);

// square periodic waveform
double squareWaveModulation(const std::array<double, 3>& pos, const __float128& t0, const __float128& t1, const __float128& beta);

#endif