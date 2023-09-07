#ifndef AISFREEPROPAGATOR_HH
#define AISFREEPROPAGATOR_HH

#include "AISAtomEnsemble.hh"
#include "AISAtom.hh"
#include "AISWavePacket.hh"

#include "AISConstants.hh"

class AISFreePropagator
{
protected:
    void PropagateAtom(AISAtom* atom, __float128 time);
    void PropagateWavePacket(AISWavePacket* wavePacket, __float128 time);
    
    virtual doubleThreeVector CalculateNewPos(doubleThreeVector pos0, doubleThreeVector vel0, __float128 time) = 0;
    virtual doubleThreeVector CalculateNewVel(doubleThreeVector pos0, doubleThreeVector vel0, __float128 time) = 0;

    doubleThreeVector newPos = {0., 0., 0.};
    doubleThreeVector newVel = {0., 0., 0.};

    doubleThreeVector currentPos = {0., 0., 0.};
    doubleThreeVector currentVel = {0., 0., 0.};
    __float128 currentTime = 0.0q;

public:
    AISFreePropagator(/* args */);
    ~AISFreePropagator();

    void PropagateEnsemble(AISAtomEnsemble* atomEnsemble, __float128 time);
};

#endif