#ifndef AISABCDPROPAGATOR_HH
#define AISABCDPROPAGATOR_HH

#include "AISFreePropagator.hh"
#include "AISUtilities.hh"

class AISAbcdPropagator : public AISFreePropagator
{
public:
    AISAbcdPropagator(__float128 dt);
    ~AISAbcdPropagator();

//protected:
    // Hamiltonian coefficients
    double3x3Matrix alpha, gamma;
    doubleThreeVector gVector;
    // ABCDXiPhi matrices 
    double3x3Matrix A, B, C, D, I;
    doubleThreeVector Xi, Phi;

    // Propagation parameters
    __float128 deltaTime;
    double deltaTime64;
    bool fAddEnergyPhase = false;
    
    // Methods to compute the ABCDXiPhi matrices
    void setABCDXiPhi();

    // Kinematic propagation
    doubleThreeVector CalculateNewPos(const doubleThreeVector& pos0, const doubleThreeVector& vel0) override;
    doubleThreeVector CalculateNewVel(const doubleThreeVector& pos0, const doubleThreeVector& vel0) override;

    // Phase propagation
    double CalculateNewPhaseDouble(const double& phase0, const doubleThreeVector& pos0, const doubleThreeVector& pos1, 
                                   const doubleThreeVector& vel0, const doubleThreeVector& vel1) override;
    
    __float128 CalculateNewPhaseQuad(const __float128& phase0) override;
};


#endif