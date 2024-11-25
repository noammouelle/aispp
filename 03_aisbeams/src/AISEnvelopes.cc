#include "AISEnvelopes.hh"

// square pulse and flat envelope
double flatSquareEnvelope(const std::array<double, 3>& pos, const __float128& t0, const __float128& t1)
{
    return 1.0;
}

// square waveform
double squareWaveModulation(const std::array<double, 3>& pos, const __float128& t0, const __float128& t1, const __float128& beta)
{
    __float128 phase = beta * (t1 - t0);
    double arg       = static_cast<double>(cosf128(phase));
    if (arg > 0) return 1;
    if (arg < 0) return -1;
    return 0;
}