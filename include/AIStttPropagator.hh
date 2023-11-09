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

class AIStttPropagator
{
private:
    __float128 fDeltaTime;
    double fDeltaTime64;
    AISIntensityProfile* fIntensityProfile;
    AISWaveFront* fWaveFront;

    __float128 fOmega;
    __float128 fOmega0;
    doubleThreeVector fK;

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
    
public:
    AIStttPropagator(__float128 deltaTime);
    ~AIStttPropagator();

    void SetWaveFront(AISWaveFront* waveFront);
    void SetIntensityProfile(AISIntensityProfile* intensityProfile);

    doubleThreeVector GetWaveVector();
    doubleThreeVector GetDoubleWaveVector();
    void SetWaveVector(const doubleThreeVector& waveVector);

    __float128 GetOmega();
    void SetOmega(const __float128& omega);

    void PropagateEnsemble(AISAtomEnsemble* atomEnsemble);
    void TttPropagateEnsemble(AISAtomEnsemble* AISAtomEnsemble);
    void TttPropagateAtom(AISAtom* atom);
    void TttPropagateWavePacket(AISWavePacket* wavePacket);
};


#endif