#include "AISEnvelopes.hh"

// square pulse and flat envelope
double flatSquareEnvelope(const std::array<double, 3>& pos, 
                          const double& param1,
                          const __float128& t0, const __float128& t1)
{
    return 1.0;
}

// gaussian envelope (waist at z=0)
double gaussianEnvelope(const std::array<double, 3>& pos, 
                        const double& w0,
                        const __float128& t0, const __float128& t1)
{
    double r   = sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
    double z   = pos[2];
    double zR  = pi * w0 * w0 / lambdaSr87;
    double w   = w0 * sqrt(1 + (z/zR)*(z/zR));

    double envelope = exp(- r * r / (w * w)) * w0 / w;

    return envelope;
}