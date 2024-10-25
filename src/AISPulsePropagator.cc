#include "AISPulsePropagator.hh"

AISPulsePropagator::AISPulsePropagator(__float128 deltaTime)
{
    fDeltaTime = deltaTime;
    fDeltaTime64 = convertScalarToDouble(fDeltaTime);
}

AISPulsePropagator::~AISPulsePropagator()
{
    // in practice those are defined externally and should not be deleted here I think
    
    //delete fIntensityProfile;
    //delete fWaveFront;
}

void AISPulsePropagator::SetWaveFront(AISWaveFront* waveFront)
{
    fWaveFront = waveFront;
}

void AISPulsePropagator::SetIntensityProfile(AISIntensityProfile* intensityProfile)
{
    fIntensityProfile = intensityProfile;
}

doubleThreeVector AISPulsePropagator::GetWaveVector()
{
    return fK;
}

void AISPulsePropagator::SetWaveVector(const doubleThreeVector& waveVector)
{
    fK = waveVector;
}

__float128 AISPulsePropagator::GetOmega()
{
    return fOmega;
}

void AISPulsePropagator::SetOmega(const __float128& omega)
{
    fOmega = omega;
}

void AISPulsePropagator::PropagateEnsemble(AISAtomEnsemble* atomEnsemble)
{
    #pragma omp parallel for
    for(int i = 0; i < atomEnsemble->GetNumberOfAtoms(); ++i){
        PropagateAtom(atomEnsemble->GetAtom(i));
    }
}

