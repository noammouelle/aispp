#ifndef AISTTTPROPAGATOR_HH
#define AISTTTPROPAGATOR_HH

#include "AISUtilities.hh"
#include "AISWaveFront.hh"
#include "AISIntensityProfile.hh"
#include "AISAtomEnsemble.hh"
#include "AISFreePropagator.hh"
#include "AISConstants.hh"
#include "AISAbcdUtilities.hh"
#include "AISLinearGravityPropagator.hh"
#include "AISAtom.hh"
#include "AISPulsePropagator.hh"

class AIStttPropagator : public AISPulsePropagator
{
protected:
    AISLinearGravityPropagator* fFreePropagator;

    // Hamiltonian coefficients (same as the one for the propagator)
    double3x3Matrix alpha, gamma;
    doubleThreeVector gVector;

    // ABCD matrices for T = + deltaTime / 2
    double3x3Matrix A, B, C, D, I; 
    doubleThreeVector Xi, Phi;
    // ABCD matrices for T = - deltaTime / 2
    double3x3Matrix A_, B_, C_, D_;
    doubleThreeVector Xi_, Phi_;

    // methods to compute detuning and rotation matrix
    double computeDetuning(const doubleThreeVector& pos, const doubleThreeVector& vel,
                           const __float128& omega);
    complexDouble Smatrix(const doubleThreeVector& pos, const doubleThreeVector& vel,
                          const __float128& omega, const double& detuning, 
                          const std::string& transition);

    // amplitude threshold
    double amplitudeThreshold;
    
public:
    AIStttPropagator(__float128 deltaTime);
    ~AIStttPropagator();
    
    void TttPropagateAtom(AISAtom* atom);

    void PropagateAtom(AISAtom* atom) override;

    void setAmplitudeThreshold(double A);
};


#endif