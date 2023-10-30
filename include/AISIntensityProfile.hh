#ifndef AISINTENSITYPROFILE_HH
#define AISINTENSITYPROFILE_HH

#include <array>
#include <cmath>
#include "AISUtilities.hh"

class AISIntensityProfile
{
private:
    /* data */
public:
    AISIntensityProfile(double BeamRadius, double centerRabiFreq);
    ~AISIntensityProfile();

    double BeamRadius;
    double centerRabiFreq;

    double GetEffectiveRabiFreq(doubleThreeVector pos);
};



#endif