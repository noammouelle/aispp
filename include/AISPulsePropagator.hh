#ifndef AISPULSEPROPAGATOR_HH
#define AISPULSEPROPAGATOR_HH

#include "AISIntensityProfile.hh"
#include "AISWaveFront.hh"
#include "AISSineWaveFront.hh"
#include "AISUtilities.hh"
#include "AISAtomEnsemble.hh"
#include "AISAtom.hh"
#include "AISWavePacket.hh"


class AISPulsePropagator
{
public:
    AISPulsePropagator(__float128 deltaTime);
    ~AISPulsePropagator();

//protected:
    /* data */
    __float128 fDeltaTime;
    double fDeltaTime64;

    AISIntensityProfile* fIntensityProfile;
    AISSineWaveFront* fWaveFront;

    __float128 fOmega;
    __float128 fOmega0;
    doubleThreeVector fK;

    /* methods */
    void SetWaveFront(AISSineWaveFront* waveFront);
    void SetIntensityProfile(AISIntensityProfile* intensityProfile);

    doubleThreeVector GetWaveVector();
    doubleThreeVector GetDoubleWaveVector();
    void SetWaveVector(const doubleThreeVector& waveVector);

    __float128 GetOmega();
    void SetOmega(const __float128& omega);

    void PropagateEnsemble(AISAtomEnsemble* atomEnsemble);

    /* virtual methods */
    virtual void PropagateAtom(AISAtom* atom) = 0;
};


#endif