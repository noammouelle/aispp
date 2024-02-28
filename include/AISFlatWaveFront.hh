#ifndef AISFLATWAVEFRONT_HH
#define AISFLATWAVEFRONT_HH

#include "AISUtilities.hh"
#include "AISWaveFront.hh"

class AISFlatWaveFront : public AISWaveFront
{
public:
    AISFlatWaveFront(doubleThreeVector psrGradient, double laserPhase);
    ~AISFlatWaveFront();

    double GetValue(const doubleThreeVector& pos) override;
};

#endif