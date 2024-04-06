#ifndef AISLMTPROPAGATOR_HH
#define AISLMTPROPAGATOR_HH

#include "AISUtilities.hh"
#include "AISAtom.hh"
#include "AISWaveFront.hh"
#include "AISSineWaveFront.hh"
#include "AISIntensityProfile.hh"

class AISLmtPropagator
{
public:
    AISLmtPropagator(int lmtOrder, double vz0, __float128 lmtDelayTime, __float128 lmtPulseTime, int lmtBlockIndex);
    ~AISLmtPropagator();

    int GetLmtOrder();
    double GetVz0();
    __float128 GetLmtDelayTime();
    int GetLmtBlockIndex();
    int GetNpulses();

    void SetWaveFronts(AISWaveFront* upwardWavefront, AISWaveFront* downwardWavefront);
    void SetIntensityProfiles(AISIntensityProfile* upwardIntensityProfile, AISIntensityProfile* downwardIntensityProfile);

    void PropagateAtom(AISAtom* atom);

    virtual double computeVelocity(const double& vz0, const __float128& t) = 0;
    virtual void PropagateAtomInPulse(AISAtom* atom, const doubleThreeVector& k, const __float128& omega,
                                      AISWaveFront* waveFront, AISIntensityProfile* intensityProfile) = 0;
    virtual void PropagateAtomFreely(AISAtom* atom) = 0;


//protected:
    int fNpulses;
    int fLmtOrder;
    __float128 fLmtDelayTime;
    int fLmtBlockIndex;
    double fVz0;

    AISWaveFront* fUpwardWavefront;
    AISWaveFront* fDownwardWavefront;
    AISIntensityProfile* fUpwardIntensityProfile;
    AISIntensityProfile* fDownwardIntensityProfile;
    doubleThreeVector fK0;

    doubleVector fAtomVelocities; // used to compute detunings
    std::vector<doubleThreeVector> fUpwardWaveVectors;
    std::vector<doubleThreeVector> fDownwardWaveVectors;
    quadVector fUpwardOmegas;
    quadVector fDownwardOmegas;

    void computeVelocities();
    void computeWaveVectors();
    void computeOmegas();
};

#endif