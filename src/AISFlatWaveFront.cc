#include "AISFlatWaveFront.hh"

AISFlatWaveFront::AISFlatWaveFront(doubleThreeVector psrGradient, double laserPhase) : AISWaveFront(psrGradient, laserPhase)
{
}

AISFlatWaveFront::~AISFlatWaveFront()
{
}

double AISFlatWaveFront::GetValue(const doubleThreeVector& pos, const doubleThreeVector& k)
{
    return dotProduct(fPsrGradient, pos) + fLaserPhase;
}