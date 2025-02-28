#ifndef AISKINEMATICPROPAGATOR_HH
#define AISKINEMATICPROPAGATOR_HH

#include "AISAtomEnsemble.hh"
#include "AISAtom.hh"
#include "AISWavePacket.hh"

#include "AISConstants.hh"

#include <stdio.h>
#include <math.h>
#include <memory>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_integration.h>

class AISKinematicPropagator
{
public: 
    using potentialFunctionType = double(*)(const doubleThreeVector&, const doubleThreeVector&);
    using gradPotentialFunctionType = doubleThreeVector(*)(const doubleThreeVector&, const doubleThreeVector&);
    using hessianPotentialFunctionType = double3x3Matrix(*)(const doubleThreeVector&, const doubleThreeVector&);

    std::string potentialTypeName;

    std::string GetPotentialTypeName();
    void SetPotentialTypeName(std::string potentialTypeName);

    AISKinematicPropagator(std::shared_ptr<potentialFunctionType> U, std::shared_ptr<gradPotentialFunctionType> dUdx,
                           std::shared_ptr<gradPotentialFunctionType> dUdp,
                           std::shared_ptr<hessianPotentialFunctionType> d2Udxdx,
                           std::shared_ptr<hessianPotentialFunctionType> d2Udxdp,
                           std::shared_ptr<hessianPotentialFunctionType> d2Udpdp);
    ~AISKinematicPropagator();

    void PropagateEnsemble(std::unique_ptr<AISAtomEnsemble>& atomEnsemble, __float128 t1);
    void PropagateAtom(std::unique_ptr<AISAtom>& atom, __float128 t1);
    void PropagateWavePacket(std::unique_ptr<AISWavePacket>& wavePacket, __float128 t1);

    void PropagateEnsembleLinearized(std::unique_ptr<AISAtomEnsemble>& atomEnsemble, __float128 t1);
    void PropagateAtomLinearized(std::unique_ptr<AISAtom>& atom, __float128 t1);
    void PropagateWavePacketLinearized(std::unique_ptr<AISWavePacket>& wavePacket, __float128 t1);

    void SetAddEnergyPhase(bool addEnergyPhase);

    std::array<double,2> CalculateNewPhaseDouble(const double& phase0, const double& phaseErr, const doubleThreeVector& pos0, const doubleThreeVector& vel0,
                                                 const doubleThreeVector& atomInitialPos, const doubleThreeVector& atomInitialVel,
                                                 const double t0, const double t1, const __float128& atomInitialTime);
    __float128 CalculateNewPhaseQuad(const __float128& phase0, const __float128& t0, const __float128 t1);

    std::array<doubleThreeVector, 2> CalculateNewPhaseSpaceCoords(const double& t0, const double& t1, 
                                                                  const doubleThreeVector& pos0, const doubleThreeVector& vel0);

    std::array<doubleThreeVector, 2> CalculateNewPhaseSpaceCoordsLinearized(const double& t0, const double& t1, 
                                                                            const doubleThreeVector& pos0, const doubleThreeVector& vel0,
                                                                            const doubleThreeVector& posStar, const doubleThreeVector& velStar); 
    std::array<doubleThreeVector, 2> get_dotPhaseSpaceCoordsLinearized(const double& t0, const double& t1, 
                                                                       const doubleThreeVector& posPrime, const doubleThreeVector& velPrime,
                                                                       const doubleThreeVector& posStar, const doubleThreeVector& velStar);

    std::tuple<double3x3Matrix, double3x3Matrix, doubleThreeVector> get_ABXi(const double& t0, const double& t1, 
                                                                             const doubleThreeVector& posStar, const doubleThreeVector& velStar);
    
    // These member functions return the Lagrangian and the action divided by m/hbar
    static double get_dL_wrapper(double t, void *params); 
    double get_dL(const double& t, void *params);
    std::array<double,2> get_dScl(const double& t0, const double& t1, const doubleThreeVector& pos0, const doubleThreeVector& vel0,
                                 const doubleThreeVector& atomInitialPos, const doubleThreeVector& atomInitialVel,
                                 const __float128& atomInitialTime);

    static int func(double t, const double y[], double f[], void *params);
    static int funcLinearized(double t, const double y[], double f[], void *params);

    doubleThreeVector get_dUdx(const doubleThreeVector& pos, const doubleThreeVector& vel);
    doubleThreeVector get_dUdp(const doubleThreeVector& pos, const doubleThreeVector& vel);
    double3x3Matrix get_d2Udxdx(const doubleThreeVector& pos, const doubleThreeVector& vel);
    double3x3Matrix get_d2Udxdp(const doubleThreeVector& pos, const doubleThreeVector& vel);
    double3x3Matrix get_d2Udpdp(const doubleThreeVector& pos, const doubleThreeVector& vel);

    // Careful! This function returns the potential energy divided by m
    double get_U(const doubleThreeVector& pos, const doubleThreeVector& vel);

    __float128 deltaTime = 0.0q;

    bool fAddEnergyPhase = false;//true;//false;

    std::shared_ptr<potentialFunctionType> U;
    std::shared_ptr<gradPotentialFunctionType> dUdx, dUdp;
    std::shared_ptr<hessianPotentialFunctionType> d2Udxdx, d2Udxdp, d2Udpdp;

    // tolerance for the integration and the ODE solvers
    double qagAbsTol = 1e-12;
    double qagRelTol = 1e-12;
    double odeAbsTol = 1e-9;
    double odeRelTol = 0.0;
};

// struct for the Lagrangian parameters
struct LagrangianParams {
    double t0;
    doubleThreeVector pos0;
    doubleThreeVector vel0;
    doubleThreeVector atomInitialPos;
    doubleThreeVector atomInitialVel;
    __float128 atomInitialTime;
    AISKinematicPropagator* propagator;
};

struct FuncLinearizedParams {
    doubleThreeVector posStar; // phase space coordinates around which the Hamiltonian is expanded
    doubleThreeVector velStar;
    AISKinematicPropagator* propagator;
    doubleThreeVector dUdx;
    doubleThreeVector dUdp;
    double3x3Matrix d2Udxdx;
    double3x3Matrix d2Udxdp;
    double3x3Matrix d2Udpdx;
    double3x3Matrix d2Udpdp;
};

#endif