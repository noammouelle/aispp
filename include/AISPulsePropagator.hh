#ifndef AISPULSEPROPAGATOR_HH
#define AISPULSEPROPAGATOR_HH

#include "AISAtomEnsemble.hh"
#include "AISAtom.hh"
#include "AISWavePacket.hh"
#include "AISLaserBeam.hh"
#include "AISKinematicPropagator.hh"

#include "AISConstants.hh"

#include <stdio.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_integration.h>

class AISPulsePropagator
{
private:
    AISLaserBeam* laserBeam;
    AISKinematicPropagator* kinematicPropagator;
    __float128 t0, t1;

protected:
    void ApplyU1(AISWavePacket* wavepacket, __float128 t0, __float128 t1);
    void ApplyU2(AISWavePacket* wavepacket, __float128 t0); // in theory could define in terms of t0 and t1 but in practice only need U2(t0,t0).
    void ApplyU2Dagger(AISWavePacket* wavepacket, __float128 t0, __float128 t1);
    void ApplyU3(AISWavePacket* wavepacket0, AISWavePacket* wavepacket1, __float128 t0, __float128 t1);

    static int funcU3(double t, const double y[], double f[], void *params);

    double getDelta(doubleThreeVector pos0, doubleThreeVector vel0, double t0, double t1, AISLaserBeam* laserBeam);
public: 
    AISPulsePropagator(AISLaserBeam* beam, __float128 t0, __float128 t1,
                       AISKinematicPropagator* kinematicPropagator);
    ~AISPulsePropagator();

    void PropagateEnsemble(AISAtomEnsemble* atomEnsemble);
    void PropagateAtom(AISAtom* atom);
};

struct U3Params
{
    double t0;
    doubleThreeVector pos, vel;
    AISLaserBeam* laserBeam;
    AISPulsePropagator* pulsePropagator;
};

#endif