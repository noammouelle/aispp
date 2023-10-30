#include <array>
#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <cstdio>

#include "AISAtomEnsemble.hh"
#include "AISFreePropagator.hh"
#include "AISLinearGravityPropagator.hh"
#include "AIStttPropagator.hh"
#include "AISUtilities.hh"
#include "AISWaveFront.hh"
#include "AISConstants.hh"
#include "AISIntensityProfile.hh"
#include "AISWaveFront.hh"

doubleThreeVector computeDetunedWaveVector(const doubleThreeVector& k0, const double& vz)
{
    doubleThreeVector detunedK = k0;
    for(int i = 0; i < 3; ++i)
    {
        if(i == 2){
            detunedK[i] = k0[i] /(1 - vz / c);
        }
    }

    return detunedK;
}

__float128 computeDetunedOmega(const __float128& omega0, const double& vz){
    return omega0 / (1 - vz / c);
}

void writeFile(std::string fName, AISAtomEnsemble* atomEnsemble)
{
    // write a file to cross check
    std::ofstream outFile(fName);

    outFile << std::fixed; // Use fixed-point notation
    outFile << std::setprecision(15);

    outFile << "State, Amplitude, X, Y, Z, VX, VY, VZ, Phase, PhaseQuad" << std::endl;
    for(int i = 0; i < atomEnsemble->GetNumberOfAtoms(); ++i)
    {
        AISAtom* currentAtom   = atomEnsemble->GetAtom(i);
        for(int j = 0; j < currentAtom->GetNumberOfWavePackets(); ++j)
        {
            AISWavePacket* currentWavePacket = currentAtom->GetWavePacket(j);
            int currentState = currentWavePacket->GetState();
            double currentAmplitude = currentWavePacket->GetAmplitude();
            doubleThreeVector currentPos = currentWavePacket->GetPosition();
            doubleThreeVector currentVel = currentWavePacket->GetVelocity();
            double currentPhaseDouble = currentWavePacket->GetPhaseDouble();
            __float128 currentPhaseQuad = currentWavePacket->GetPhaseQuad();

            long double currentPhaseQuadAsLongDouble = static_cast<long double>(currentPhaseQuad);

            char currentPhaseQuadStr[50];  // Adjust the buffer size as needed
            snprintf(currentPhaseQuadStr, sizeof(currentPhaseQuadStr), "%.34Le", currentPhaseQuadAsLongDouble);

            outFile << currentState << "," << "," << currentAmplitude << currentPos[0] << "," << currentPos[1] << "," << currentPos[2] << "," 
                    << currentVel[0] << "," << currentVel[1] << "," << currentVel[2] << ","
                    << currentPhaseDouble << "," << currentPhaseQuadStr << std::endl;
        }
    }

    outFile.close();
}

int main()
{   
    /* CLOUD PARAMETERS */
    int numAtoms      = 100;
    double cloudTemp  = 15e-9 *0;
    double cloudWidth = 1e-3 *0;

    doubleThreeVector initialPos = {0.,0.,0.};
    doubleThreeVector initialVel = {0.,0.,0.};

    /* SEQUENCE PARAMETERS */
    __float128 dt = 0.05q;
    int nsteps = 20;

    /* PULSE PARAMETERS */
    doubleThreeVector k0 = {0., 0., omegaSr87 / c};
    doubleThreeVector detunedK = computeDetunedWaveVector(k0, - g * nsteps * convertScalarToDouble(dt));
    __float128 detunedFrequency = computeDetunedOmega(omegaSr87, - g * nsteps * dt);
    
    double beamRadius = 0.5e-2;
    double rabiFreq   = 2 * pi * 3e6;

    __float128 pulseDuration = pi128 / (2 * rabiFreq);

    doubleThreeVector psrGradient = {0., 0., 0.};
    double laserPhase  = 0.0;

    double (*aberrationFunction)(const doubleThreeVector&) = zeroAberrationFunction;
    double (*aberrationFunctionDerivative)(const doubleThreeVector&) = zeroAberrationFunction;

    /* Define all the A.I. objects */
    AISAtomEnsemble* atomEnsemble          = new AISAtomEnsemble(numAtoms, cloudTemp, cloudWidth, initialPos, initialVel);

    AISLinearGravityPropagator* propagator = new AISLinearGravityPropagator(dt);

    AISIntensityProfile* intensityProfile  = new AISIntensityProfile(beamRadius, rabiFreq);
    AISWaveFront* waveFront                = new AISWaveFront(aberrationFunction, aberrationFunctionDerivative, psrGradient, laserPhase);

    AIStttPropagator* tttPropagator        = new AIStttPropagator(pulseDuration);
    tttPropagator->SetWaveFront(waveFront);
    tttPropagator->SetIntensityProfile(intensityProfile);
    tttPropagator->SetWaveVector(detunedK);
    tttPropagator->SetOmega(detunedFrequency);

    //writeFile("/home/noammouelle/sim/data_ais++/output0.csv",atomEnsemble);
    std::cout << "Initial state written to file" << std::endl;

    // free propagation
    for(int i = 0; i < nsteps; ++i){
        propagator->PropagateEnsemble(atomEnsemble,dt);
    }

    std::cout << "Free propagation done" << std::endl;

    //writeFile("/home/noammouelle/sim/data_ais++/output1.csv",atomEnsemble);

    std::cout << "Post-free propagation state written to file" << std::endl;

    // beam-splitting
    tttPropagator->PropagateEnsemble(atomEnsemble);
    std::cout << "Beam-splitting done" << std::endl;

    //writeFile("/home/noammouelle/sim/data_ais++/output2.csv",atomEnsemble);
    std::cout << "Post-beam-splitting state written to file" << std::endl;
    
    return 0;
}






