#ifndef AISABCDPROPAGATOR_HH
#define AISABCDPROPAGATOR_HH

#include "AISFreePropagator.hh"
#include "AISUtilities.hh"

class AISAbcdPropagator : public AISFreePropagator
{
protected:
    // Hamiltonian coefficients
    quad3x3Matrix alpha, gamma;
    quadThreeVector gVector;
    // ABCDXiPhi matrices 
    quad3x3Matrix A, B, C, D, I;
    quadThreeVector Xi, Phi;
    double3x3Matrix A64, B64, C64, D64;
    doubleThreeVector Xi64, Phi64;

    // Propagation parameters
    __float128 deltaTime;
    double deltaTime64;
    
    // Methods to compute the ABCDXiPhi matrices
    void setABCDXiPhi();

    // Kinematic propagation
    doubleThreeVector CalculateNewPos(doubleThreeVector pos0, doubleThreeVector vel0, __float128 time) override;
    doubleThreeVector CalculateNewVel(doubleThreeVector pos0, doubleThreeVector vel0, __float128 time) override;

public:
    AISAbcdPropagator(__float128 dt);
    ~AISAbcdPropagator();
};


#endif