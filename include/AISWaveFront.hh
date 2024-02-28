#ifndef AISWAVEFRONT_HH
#define AISWAVEFRONT_HH

#include <array>
#include "AISUtilities.hh"

class AISWaveFront
{
public:
    doubleThreeVector fPsrGradient;
    double fLaserPhase;

    AISWaveFront(doubleThreeVector psrGradient, double laserPhase);
    ~AISWaveFront();

    virtual double GetValue(const doubleThreeVector& pos) = 0;
};


#endif