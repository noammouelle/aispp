#include "AISSineWaveFront.hh"

AISSineWaveFront::AISSineWaveFront(doubleThreeVector psrGradient, double laserPhase, double A, doubleThreeVector k) : AISWaveFront(psrGradient, laserPhase)
{
    fPsrGradient = psrGradient;
    fLaserPhase = laserPhase;
    faberrationAmplitude = A;
    faberrationK = k; // not the laser wavevector !!!
}

AISSineWaveFront::~AISSineWaveFront(){}

double AISSineWaveFront::GetValue(const doubleThreeVector& pos, const doubleThreeVector& k) 
{
    double linearPart = dotProduct(fPsrGradient, pos) + fLaserPhase;
    double aberration =  faberrationAmplitude * sin(dotProduct(faberrationK, pos));

    return linearPart + aberration;
}