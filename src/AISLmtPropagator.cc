#include "AISLmtPropagator.hh"

AISLmtPropagator::AISLmtPropagator(int lmtOrder, double vz0, __float128 lmtDelayTime, __float128 lmtPulseTime, int lmtBlockIndex)
{
    // Check that the lmt order is odd
    if(lmtOrder%2 == 0)
    {
        std::cout << "ERROR: lmt order must be odd" << std::endl;
        exit(1);
    }
    // Initialize member variables
    fLmtOrder = lmtOrder;
    fVz0 = vz0;
    fLmtDelayTime = lmtDelayTime;
    fLmtBlockIndex = lmtBlockIndex;
    fNpulses = static_cast<int>(lmtOrder/2.0 - 0.5);
}

AISLmtPropagator::~AISLmtPropagator()
{
}

int AISLmtPropagator::GetLmtOrder()
{
    return fLmtOrder;
}

double AISLmtPropagator::GetVz0()
{
    return fVz0;
}

__float128 AISLmtPropagator::GetLmtDelayTime()
{
    return fLmtDelayTime;
}

int AISLmtPropagator::GetLmtBlockIndex()
{
    return fLmtBlockIndex;
}

int AISLmtPropagator::GetNpulses()
{
    return fNpulses;
}

void AISLmtPropagator::SetWaveFronts(AISWaveFront* upwardWavefront, AISWaveFront* downwardWavefront)
{
    fUpwardWavefront = upwardWavefront;
    fDownwardWavefront = downwardWavefront;
}

void AISLmtPropagator::SetIntensityProfiles(AISIntensityProfile* upwardIntensityProfile, AISIntensityProfile* downwardIntensityProfile)
{
    fUpwardIntensityProfile = upwardIntensityProfile;
    fDownwardIntensityProfile = downwardIntensityProfile;
}

void AISLmtPropagator::computeVelocities()
{
    for(int pulseIndex = 0; pulseIndex < fNpulses; ++pulseIndex)
    {
        double currentVelocity = computeVelocity(fVz0, fLmtDelayTime * pulseIndex);
        fAtomVelocities.push_back(currentVelocity);
    }
}

void AISLmtPropagator::computeWaveVectors()
{
    if(fLmtBlockIndex % 2 == 0)
    {
        // blocks with even index always start with a downward pulse
        for(int pulseIndex = 0; pulseIndex < fNpulses; ++pulseIndex)
        {
            if(pulseIndex%2 == 0){
                fDownwardWaveVectors.push_back(computeDetunedWaveVector(fK0, - fAtomVelocities[pulseIndex])); // note the minus sign
                fDownwardOmegas.push_back(computeDetunedOmega(omegaSr87, - fAtomVelocities[pulseIndex]));     // (z axis is pointing up)
            }
            else{
                fUpwardWaveVectors.push_back(computeDetunedWaveVector(fK0, fAtomVelocities[pulseIndex]));
                fUpwardOmegas.push_back(computeDetunedOmega(omegaSr87, fAtomVelocities[pulseIndex]));
            }
        }
    } else if(fLmtBlockIndex % 2 == 1){
    // blocks with odd index start with a downward pulse if the number of pulses is odd
    // and with an upward pulse if the number of pulses is even
    if(fNpulses%2 == 0){
        for(int pulseIndex = 0; pulseIndex < fNpulses; ++pulseIndex)
        {
            if(pulseIndex%2 == 0){
                fUpwardWaveVectors.push_back(computeDetunedWaveVector(fK0, fAtomVelocities[pulseIndex]));
                fUpwardOmegas.push_back(computeDetunedOmega(omegaSr87, fAtomVelocities[pulseIndex]));
            }
            else{
                fDownwardWaveVectors.push_back(computeDetunedWaveVector(fK0, - fAtomVelocities[pulseIndex])); // note the minus sign
                fDownwardOmegas.push_back(computeDetunedOmega(omegaSr87, - fAtomVelocities[pulseIndex]));     // (z axis is pointing up)
            }
        }
    } else if(fNpulses%2 == 1){
        for(int pulseIndex = 0; pulseIndex < fNpulses; ++pulseIndex)
        {
            if(pulseIndex%2 == 0){
                fDownwardWaveVectors.push_back(computeDetunedWaveVector(fK0, - fAtomVelocities[pulseIndex])); // note the minus sign
                fDownwardOmegas.push_back(computeDetunedOmega(omegaSr87, - fAtomVelocities[pulseIndex]));     // (z axis is pointing up)
            }
            else{
                fUpwardWaveVectors.push_back(computeDetunedWaveVector(fK0, fAtomVelocities[pulseIndex]));
                fUpwardOmegas.push_back(computeDetunedOmega(omegaSr87, fAtomVelocities[pulseIndex]));
            }
        }
    }
    }
}

void AISLmtPropagator::computeOmegas()
{
    if(fLmtBlockIndex % 2 == 0){
        // blocks with even index always start with a downward pulse
        for(int pulseIndex = 0; pulseIndex < fNpulses; ++pulseIndex)
        {
            if(pulseIndex%2 == 0){
                fDownwardOmegas.push_back(computeDetunedOmega(omegaSr87, - fAtomVelocities[pulseIndex])); // (z axis is pointing up)
            }
            else{
                fUpwardOmegas.push_back(computeDetunedOmega(omegaSr87, fAtomVelocities[pulseIndex]));
            }
        
        }
    }
    else if(fLmtBlockIndex % 2 == 1){
    // blocks with odd index start with a downward pulse if the number of pulses is odd
    // and with an upward pulse if the number of pulses is even
    if(fNpulses%2 == 0){
        for(int pulseIndex = 0; pulseIndex < fNpulses; ++pulseIndex)
        {
            if(pulseIndex%2 == 0){
                fUpwardOmegas.push_back(computeDetunedOmega(omegaSr87, fAtomVelocities[pulseIndex]));
            }
            else{
                fDownwardOmegas.push_back(computeDetunedOmega(omegaSr87, - fAtomVelocities[pulseIndex])); // (z axis is pointing up)
            }
        }
    } else if(fNpulses%2 == 1){
        for(int pulseIndex = 0; pulseIndex < fNpulses; ++pulseIndex)
        {
            if(pulseIndex%2 == 0){
                fDownwardOmegas.push_back(computeDetunedOmega(omegaSr87, - fAtomVelocities[pulseIndex])); // (z axis is pointing up)
            }
            else{
                fUpwardOmegas.push_back(computeDetunedOmega(omegaSr87, fAtomVelocities[pulseIndex]));
            }
        }
    }
    }
}

void AISLmtPropagator::PropagateAtom(AISAtom* atom)
{
    if(fLmtBlockIndex % 2 == 0)
    {
        // blocks with even index always start with a downward pulse
        for(int pulseIndex = 0; pulseIndex < fNpulses; ++pulseIndex)
        {   
            // free propagation between pulses
            PropagateAtomFreely(atom);
            // propagate in pulse
            if(pulseIndex%2 == 0){
                int index = static_cast<int>(pulseIndex/2.0);
                PropagateAtomInPulse(atom, fDownwardWaveVectors[index], fDownwardOmegas[index],
                                    fDownwardWavefront, fDownwardIntensityProfile);
            }
            else if(pulseIndex%2 == 1){
                int index = static_cast<int>((pulseIndex - 1)/2.0);
                PropagateAtomInPulse(atom, fUpwardWaveVectors[index], fUpwardOmegas[index],
                                     fUpwardWavefront, fUpwardIntensityProfile);
            }
        }
    } else if(fLmtBlockIndex % 2 == 1){
    // blocks with odd index start with a downward pulse if the number of pulses is odd
    // and with an upward pulse if the number of pulses is even
    if(fNpulses%2 == 0){
        for(int pulseIndex = 0; pulseIndex < fNpulses; ++pulseIndex)
        {
            // free propagation between pulses
            PropagateAtomFreely(atom);
            // propagate in pulse
            if(pulseIndex%2 == 0){
                int index = static_cast<int>(pulseIndex/2.0);
                PropagateAtomInPulse(atom, fUpwardWaveVectors[index], fUpwardOmegas[index],
                                     fUpwardWavefront, fUpwardIntensityProfile);
            }
            else if(pulseIndex%2 == 1){
                int index = static_cast<int>((pulseIndex - 1)/2.0);
                PropagateAtomInPulse(atom, fDownwardWaveVectors[index], fDownwardOmegas[index],
                                     fDownwardWavefront, fDownwardIntensityProfile);
            }
        }
    } else if(fNpulses%2 == 1){
        for(int pulseIndex = 0; pulseIndex < fNpulses; ++pulseIndex)
        {
            // free propagation between pulses
            PropagateAtomFreely(atom);
            // propagate in pulse
            if(pulseIndex%2 == 0){
                int index = static_cast<int>(pulseIndex/2.0);
                PropagateAtomInPulse(atom, fDownwardWaveVectors[index], fDownwardOmegas[index],
                                     fDownwardWavefront, fDownwardIntensityProfile);
            }
            else if(pulseIndex%2 == 1){
                int index = static_cast<int>((pulseIndex - 1)/2.0);
                PropagateAtomInPulse(atom, fUpwardWaveVectors[index], fUpwardOmegas[index],
                                     fUpwardWavefront, fUpwardIntensityProfile);
            }
        }
    }
    }
          
}