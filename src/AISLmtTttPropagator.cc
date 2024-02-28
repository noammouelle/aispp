#include "AISLmtTttPropagator.hh"

AISLmtTttPropagator::AISLmtTttPropagator(int lmtOrder, double vz0, __float128 lmtDelayTime, __float128 lmtPulseTime, int lmtBlockIndex) 
                                        : AISLmtPropagator(lmtOrder, vz0, lmtDelayTime, lmtPulseTime, lmtBlockIndex), fLinearGravityPropagator(new AISLinearGravityPropagator(lmtDelayTime)), fTttPropagator(new AIStttPropagator(lmtPulseTime))
{
    fLinearGravityPropagator->SetAddEnergyPhase(true); 

    // Compute the velocities of the atoms at the beginning of the pulse
    computeVelocities();
    // Compute the corresponding detuned wavevectors
    computeWaveVectors();
}

AISLmtTttPropagator::~AISLmtTttPropagator()
{
    delete fLinearGravityPropagator;
    delete fTttPropagator;
}

double AISLmtTttPropagator::computeVelocity(const double& vz0, const __float128& t) 
{
    return vz0 - g * t;
}

void AISLmtTttPropagator::PropagateAtomInPulse(AISAtom* atom, const doubleThreeVector& k, const __float128& omega,
                                               AISSineWaveFront* waveFront, AISIntensityProfile* intensityProfile)
{
    fTttPropagator->SetWaveVector(k);
    fTttPropagator->SetOmega(omega);
    fTttPropagator->SetWaveFront(waveFront);
    fTttPropagator->SetIntensityProfile(intensityProfile);

    fTttPropagator->PropagateAtom(atom);
}

void AISLmtTttPropagator::PropagateAtomFreely(AISAtom* atom)
{
    fLinearGravityPropagator->PropagateAtom(atom);
}

void AISLmtTttPropagator::setAmplitudeThreshold(double A)
{
    //fTttPropagator->setAmplitudeThreshold(A);
}