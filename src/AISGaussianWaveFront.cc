#include "AISGaussianWaveFront.hh"

AISGaussianWaveFront::AISGaussianWaveFront(doubleThreeVector psrGradient, double laserPhase, double beamRadius) : AISWaveFront(psrGradient, laserPhase)
{
    fPsrGradient = psrGradient;
    fLaserPhase = laserPhase;
    this->beamRadius = beamRadius;
    w0 = beamRadius / 2; // beam waist (only param). Defined such that it is half the beam radius (ARBITRARY)
    zR = pi * pow(w0, 2) / lambdaSr87;
}

AISGaussianWaveFront::~AISGaussianWaveFront()
{
}

double AISGaussianWaveFront::R(const double& z)
{
    return z * (1 + pow(zR, 2) / pow(z, 2));
}

double AISGaussianWaveFront::w(const double& z)
{
    return w0 * sqrt(1 + pow(z, 2) / pow(zR, 2));
}

double AISGaussianWaveFront::gouyPhase(const double& z)
{
    return atan(z / zR);
}

double AISGaussianWaveFront::GetNonShiftedValue(const doubleThreeVector& pos, const doubleThreeVector& k)
{
    // get the coordinates, accounting for the tilt
    // WARNING: this assumes only rotations about the y-axis
    double tiltAngle = fPsrGradient[0] / k[2];
    double x = pos[0] - tiltAngle * pos[2];
    double y = pos[1];
    double z = pos[2] + tiltAngle * pos[0];

    double r = sqrt(pow(x, 2) + pow(y, 2));
    
    return - fPsrGradient[0] * pos[0] - k[2] * pow(r, 2) / (2 * R(z)) + gouyPhase(z) + fLaserPhase;
}

double AISGaussianWaveFront::GetValue(const doubleThreeVector& pos, const doubleThreeVector& k)
{
    return GetNonShiftedValue({pos[0], pos[1], pos[2] + zR*sqrt(3)}, k);
}