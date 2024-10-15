#ifndef AISPULSEPROPAGATOR_HH
#define AISPULSEPROPAGATOR_HH

#include "AISAtomEnsemble.hh"
#include "AISAtom.hh"
#include "AISWavePacket.hh"

#include "AISConstants.hh"

#include <stdio.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_integration.h>

class AISPulsePropagator
{
protected:
    std::array<double,2> CalculateNewPhaseDouble(const double& phase0, const double& phaseErr, const doubleThreeVector& pos0, const doubleThreeVector& vel0,
                                                         const double t0, const double t1);
    __float128 CalculateNewPhaseQuad(const __float128& phase0, const __float128& t0, const __float128 t1);

    std::array<doubleThreeVector, 2> CalculateNewPhaseSpaceCoords(const double& t0, const double& t1, 
                                                                  const doubleThreeVector& pos0, const doubleThreeVector& pos1, 
                                                                  const doubleThreeVector& vel0, const doubleThreeVector& vel1); 
    
    // These member functions return the Lagrangian and the action divided by m/hbar
    static double get_L_wrapper(double t, void *params); 
    double get_L(const double& t, void *params);
    std::array<double,2> get_Scl(const double& t0, const double& t1, const doubleThreeVector& pos0, const doubleThreeVector& vel0);

    static int func(double t, const double y[], double f[], void *params);

    virtual doubleThreeVector get_dUdx(const doubleThreeVector& pos, const doubleThreeVector& vel) = 0;
    virtual doubleThreeVector get_dUdp(const doubleThreeVector& pos, const doubleThreeVector& vel) = 0;

    // This function returns the potential energy divided by m
    virtual double get_U(const doubleThreeVector& pos, const doubleThreeVector& vel) = 0;

    __float128 deltaTime = 0.0q;

    bool fAddEnergyPhase = false;

    void ApplyU1(AISWavePacket* wavepacket, __float128 t0, __float128 t1);
    void ApplyU2(AISWavePacket* wavepacket, __float128 t0, __float128 t1);

public: 
    AISPulsePropagator(doubleThreeVector k, __float128 t0, __float128 t1);
    ~AISPulsePropagator();

    void PropagateEnsemble(AISAtomEnsemble* atomEnsemble, __float128 t1);
    void PropagateAtom(AISAtom* atom, __float128 t1);
    void PropagateWavePacket(AISWavePacket* wavePacket, __float128 t1);
};

#endif