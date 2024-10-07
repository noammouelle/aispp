#ifndef AISKINEMATICPROPAGATOR_HH
#define AISKINEMATICPROPAGATOR_HH

#include "AISAtomEnsemble.hh"
#include "AISAtom.hh"
#include "AISWavePacket.hh"

#include "AISConstants.hh"

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <math.h>
#include <gsl/gsl_integration.h>

class AISKinematicPropagator
{
protected:
    std::array<double,2> CalculateNewPhaseDouble(const double& phase0, const doubleThreeVector& pos0, const doubleThreeVector& vel0,
                                                         const double t0, const double t1);
    virtual __float128 CalculateNewPhaseQuad(const __float128& phase0, const __float128& t0, const __float128 t1) = 0;

    std::array<doubleThreeVector, 2> CalculateNewPhaseSpaceCoords(const double& t0, const double& t1, 
                                                                  const doubleThreeVector& pos0, const doubleThreeVector& pos1, 
                                                                  const doubleThreeVector& vel0, const doubleThreeVector& vel1) = 0;  
    double get_L(const double& t, void *params);
    std::array<double,2> get_Scl(const double& t0, const double& t1, const doubleThreeVector& pos0, const doubleThreeVector& vel0);
    int func(double t, const double y[], double f[], void *params);

    virtual doubleThreeVector get_dUdx(const doubleThreeVector& pos, const doubleThreeVector& vel) = 0;
    virtual doubleThreeVector get_dUdp(const doubleThreeVector& pos, const doubleThreeVector& vel) = 0;
    virtual double get_U(const doubleThreeVector& pos, const doubleThreeVector& vel) = 0;

    __float128 deltaTime = 0.0q;

    bool fAddEnergyPhase = false;

public:
    AISKinematicPropagator();
    ~AISKinematicPropagator();

    void PropagateEnsemble(AISAtomEnsemble* atomEnsemble, __float128 t1);
    void PropagateAtom(AISAtom* atom, __float128 t1);
    void PropagateWavePacket(AISWavePacket* wavePacket, __float128 t1);

    void SetAddEnergyPhase(bool addEnergyPhase);
};

// struct for the Lagrangian parameters
struct LagrangianParams {
    double t0;
    doubleThreeVector pos0;
    doubleThreeVector vel0;
};

#endif