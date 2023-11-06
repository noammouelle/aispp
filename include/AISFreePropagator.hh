#ifndef AISFREEPROPAGATOR_HH
#define AISFREEPROPAGATOR_HH

#include "AISAtomEnsemble.hh"
#include "AISAtom.hh"
#include "AISWavePacket.hh"

#include "AISConstants.hh"

class AISFreePropagator
{
protected:
    void PropagateAtom(AISAtom* atom);
    void PropagateWavePacket(AISWavePacket* wavePacket);
    
    virtual doubleThreeVector CalculateNewPos(const doubleThreeVector& pos0, const doubleThreeVector& vel0) = 0;
    virtual doubleThreeVector CalculateNewVel(const doubleThreeVector& pos0, const doubleThreeVector& vel0) = 0;
    virtual double CalculateNewPhaseDouble(const double& phase0, const doubleThreeVector& pos0, const doubleThreeVector& pos1, 
                                           const doubleThreeVector& vel0, const doubleThreeVector& vel1) = 0;
    virtual __float128 CalculateNewPhaseQuad(const __float128& phase0) = 0;                                      

    doubleThreeVector newPos = {0., 0., 0.};
    doubleThreeVector newVel = {0., 0., 0.};

    doubleThreeVector currentPos = {0., 0., 0.};
    doubleThreeVector currentVel = {0., 0., 0.};

    double currentPhaseDouble   = 0.0;
    __float128 currentPhaseQuad = 0.0q;

    __float128 currentTime = 0.0q;
    __float128 deltaTime = 0.0q;

    bool fAddEnergyPhase = false;

public:
    AISFreePropagator(__float128 dt);
    ~AISFreePropagator();

    void PropagateEnsemble(AISAtomEnsemble* atomEnsemble);

    void SetAddEnergyPhase(bool addEnergyPhase);
};

#endif