#include "AISFlatWaveFront.hh"

AISFlatWaveFront::AISFlatWaveFront(doubleThreeVector psrGradient, double laserPhase) : AISWaveFront(psrGradient, laserPhase)
{
}

AISFlatWaveFront::~AISFlatWaveFront()
{
}

double AISFlatWaveFront::GetValue(const doubleThreeVector& pos)
{
    return dotProduct(fPsrGradient, pos) + fLaserPhase;
}