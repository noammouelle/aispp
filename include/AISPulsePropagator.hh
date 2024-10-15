#ifndef AISPULSEPROPAGATOR_HH
#define AISPULSEPROPAGATOR_HH

#include "AISAtomEnsemble.hh"
#include "AISAtom.hh"
#include "AISWavePacket.hh"
#include "AISLaserBeam.hh"

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
    __float128 t0, t1;

protected:
    void ApplyU1(AISWavePacket* wavepacket, __float128 t0, __float128 t1);
    void ApplyU2(AISWavePacket* wavepacket, __float128 t0, __float128 t1);
    void ApplyU2Dagger(AISWavePacket* wavepacket, __float128 t0, __float128 t1);
    void ApplyU3(AISWavePacket* wavepacket0, AISWavePacket* wavepacket1, __float128 t0, __float128 t1);

public: 
    AISPulsePropagator(AISLaserBeam* beam, __float128 t0, __float128 t1);
    ~AISPulsePropagator();

    void PropagateEnsemble(AISAtomEnsemble* atomEnsemble);
    void PropagateAtom(AISAtom* atom);
};

#endif