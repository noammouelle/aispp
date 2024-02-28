#ifndef AISLMTTTTPROPAGATOR_HH
#define AISLMTTTTPROPAGATOR_HH

#include "AISLmtPropagator.hh"	
#include "AISWaveFront.hh"
#include "AISIntensityProfile.hh"
#include "AISAtom.hh"
#include "AISUtilities.hh"
#include "AISLinearGravityPropagator.hh"
#include "AIStttPropagator.hh"

class AISLmtTttPropagator : public AISLmtPropagator
{
public:
    AISLmtTttPropagator(int lmtOrder, double vz0, __float128 lmtDelayTime, __float128 lmtPulseTime, int lmtBlockIndex);
    ~AISLmtTttPropagator();

    double computeVelocity(const double& vz0, const __float128& t) override;
    void PropagateAtomInPulse(AISAtom* atom, const doubleThreeVector& k, const __float128& omega,
                              AISSineWaveFront* waveFront, AISIntensityProfile* intensityProfile) override;
    void PropagateAtomFreely(AISAtom* atom) override;

    void setAmplitudeThreshold(double A);

protected:
    AISLinearGravityPropagator* fLinearGravityPropagator;
    AIStttPropagator* fTttPropagator;
};

#endif