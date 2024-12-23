#include "AISWavefronts.hh"

// flat wavefront and gradient wavefront
double flatWavefront(const std::array<double, 3>& pos, const double& param1, const double& param2)
{
    return 0.0;
}

std::array<double, 3> flatGradientWavefront(const std::array<double, 3>& pos, const double& paaram1, const double& param2)
{
    return {0.0, 0.0, 0.0};
}

// gaussian wavefront and gradient
double gaussianWavefront(const std::array<double, 3>& pos, const double& kz,const double& w0)
{
    double r   = sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
    double z   = pos[2];
    double zR  = pi * w0 * w0 / lambdaSr87;
    double RInv = z / (z*z + zR*zR);

    return kz * r * r * RInv / 2.0 - atan(z/zR);
}

std::array<double,3> gaussianGradientWavefront(const std::array<double, 3>& pos, const double& kz, const double& w0)
{
    double r   = sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
    double x   = pos[0];
    double y   = pos[1];
    double z   = pos[2];
    double zR  = pi * w0 * w0 / lambdaSr87;
    double RInv = z / (z*z + zR*zR);

    double dWdX = kz * x * RInv;
    double dWdY = kz * y * RInv;
    double dWdZ = r * r * RInv * RInv / 2.0 - zR / (z*z + zR*zR);

    return {dWdX, dWdY, dWdZ};
}