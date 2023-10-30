#include "AISIntensityProfile.hh"

AISIntensityProfile::AISIntensityProfile(double BeamRadius, double centerRabiFreq)
{
    this->BeamRadius = BeamRadius;
    this->centerRabiFreq = centerRabiFreq;
}

AISIntensityProfile::~AISIntensityProfile()
{
}

double AISIntensityProfile::GetEffectiveRabiFreq(doubleThreeVector pos)
{
    double r = sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
    double effectiveRabiFreq = centerRabiFreq * exp(-2 * r * r / (BeamRadius * BeamRadius));
    return effectiveRabiFreq;
}