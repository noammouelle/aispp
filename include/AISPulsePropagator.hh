#ifndef AISPULSEPROPAGATOR_HH
#define AISPULSEPROPAGATOR_HH

#include "AISIntensityProfile.hh"
#include "AISWaveFront.hh"
#include "AISUtilities.hh"


class AISPulsePropagator
{
public:
    AISPulsePropagator();
    ~AISPulsePropagator();

protected:
    __float128 fDeltaTime;
    double fDeltaTime64;

    AISIntensityProfile* fIntensityProfile;
    AISWaveFront* fWaveFront;

    __float128 fOmega;
    __float128 fOmega0;
    doubleThreeVector fK;
};


#endif