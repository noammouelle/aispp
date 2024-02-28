#include "AISSineWaveFront.hh"

AISSineWaveFront::AISSineWaveFront(doubleThreeVector psrGradient, double laserPhase, double A, doubleThreeVector k) : AISWaveFront(psrGradient, laserPhase)
{
    fPsrGradient = psrGradient;
    fLaserPhase = laserPhase;
    faberrationAmplitude = A;
    faberrationK = k;
}

AISSineWaveFront::~AISSineWaveFront(){}

double AISSineWaveFront::GetValue(const doubleThreeVector& pos) 
{
    double linearPart = dotProduct(fPsrGradient, pos) + fLaserPhase;
    double aberration =  faberrationAmplitude * sin(dotProduct(faberrationK, pos));

    return linearPart + aberration;
}