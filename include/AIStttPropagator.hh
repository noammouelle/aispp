#ifndef AISTTTPROPAGATOR_HH
#define AISTTTPROPAGATOR_HH

#include <array>

using threeQuadVector = std::array<__float128,3>;

#include "AISWaveFront.hh"
#include "AISInstensityProfile.hh"
#include "AISAtomEnsemble.hh"
#include "AISFreePropagator.hh"
#include "AISConstants.hh"

class AIStttPropagator
{
private:
    __float128 fDeltaTime;
    AISInstensityProfile* fIntensityProfile;
    AISWaveFront* fWaveFront;

    __float128 fOmega;
    threeQuadVector fK;
    threeVector fDoubleK;

    AISFreePropagator* fFreePropagator;

    // some dummy variables
    threeVector dx = {0., 0., 0.};
    threeVector dv = {0., 0., 0.};
    
public:
    AIStttPropagator(__float128 deltaTime);
    ~AIStttPropagator();

    void SetWaveFront(AISWaveFront* waveFront);
    void SetIntensityProfile(AISInstensityProfile* intensityProfile);

    threeQuadVector GetWaveVector();
    threeVector GetDoubleWaveVector();
    void SetWaveVector(threeQuadVector waveVector);

    __float128 GetOmega();
    void SetOmega(__float128 omega);

    void PropagateEnsemble(AISAtomEnsemble* atomEnsemble);
    void TttPropagateEnsemble(AISAtomEnsemble* AISAtomEnsemble);
    void TttPropagateAtom(AISAtom* atom);
    void TttPropagateWavePacket(AISWavePacket* wavePacket);
};


#endif