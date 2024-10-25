#ifndef AISSINEWAVEFRONT_HH
#define AISSINEWAVEFRONT_HH

#include "AISWaveFront.hh"
#include "AISUtilities.hh"

class AISSineWaveFront : public AISWaveFront
{
public:
    double faberrationAmplitude;
    doubleThreeVector faberrationK;

    AISSineWaveFront(doubleThreeVector psrGradient, double laserPhase, double A, doubleThreeVector k);
    ~AISSineWaveFront();

    double GetValue(const doubleThreeVector& pos, const doubleThreeVector& k) override;
};

#endif