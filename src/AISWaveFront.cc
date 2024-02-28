#include "AISWaveFront.hh"
#include "AISUtilities.hh"
#include "AISConstants.hh"

// PARENT CLASS
AISWaveFront::AISWaveFront(doubleThreeVector psrGradient, double laserPhase)
{
    fPsrGradient = psrGradient;
    fLaserPhase = laserPhase;
}

AISWaveFront::~AISWaveFront(){}
