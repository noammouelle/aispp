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
#include <memory>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_integration.h>

class AISPulsePropagator
{
//private:
public:
    std::shared_ptr<AISLaserBeam> laserBeam;
    std::shared_ptr<AISKinematicPropagator> kinematicPropagator;
    __float128 initTime, finalTime;
    double initTimeDouble, finalTimeDouble;

protected:
    void ApplyU1(std::unique_ptr<AISWavePacket>& wavepacket, __float128 t0, __float128 t1);
    void ApplyU2(std::unique_ptr<AISWavePacket>& wavepacket, __float128 t0); // in theory could define in terms of t0 and t1 but in practice only need U2(t0,t0).
    void ApplyU2Dagger(std::unique_ptr<AISWavePacket>& wavepacket, __float128 t0, __float128 t1);
    void ApplyU3(std::unique_ptr<AISWavePacket>& wavepacket0, std::unique_ptr<AISWavePacket>& wavepacket1, __float128 t0, __float128 t1);

    static int funcU3(double t, const double y[], double f[], void *params);

    double getDelta(doubleThreeVector pos0, doubleThreeVector vel0, double t0, double t1, std::shared_ptr<AISLaserBeam> laserBeam);
public: 
    AISPulsePropagator(std::shared_ptr<AISLaserBeam> beam, __float128 t0, __float128 t1,
                       std::shared_ptr<AISKinematicPropagator> kinematicPropagator);
    ~AISPulsePropagator();

    void PropagateEnsemble(std::unique_ptr<AISAtomEnsemble>& atomEnsemble);
    void PropagateAtom(std::unique_ptr<AISAtom>& atom);
};

struct U3Params
{
    double t0;
    doubleThreeVector pos, vel;
    std::shared_ptr<AISLaserBeam> laserBeam;
    AISPulsePropagator* pulsePropagator;
};

#endif