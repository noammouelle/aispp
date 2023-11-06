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
#include "AISAtomInterferometer.hh"
#include "AISAtomInterferometerParams.hh"

void writeFile(std::string fName, AISAtomEnsemble* atomEnsemble)
{
    // write a file to cross check
    std::ofstream outFile(fName);
    std::cout << "Writing to file: " << fName << std::endl;
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

            outFile << currentState << "," << currentAmplitude << "," << currentPos[0] << "," << currentPos[1] << "," << currentPos[2] << "," 
                    << currentVel[0] << "," << currentVel[1] << "," << currentVel[2] << ","
                    << currentPhaseDouble << "," << currentPhaseQuadStr << std::endl;
        }
    }

    outFile.close();
}

int main()
{   
    // create a parameter object
    AISAtomInterferometerParams params;

    /* CLOUD PARAMETERS */
    params.nAtoms          = 1;
    params.cloudTemperature  = 15e-9 *0;
    params.cloudWidth        = 1e-3 *0;

    params.initialPosition   = {0.,0.,0.};
    params.initialVelocity   = {0.,0.,0.};

    /* SEQUENCE PARAMETERS */
    params.interogationTime  = 1.0q;
    params.nSteps            = 20;

    /* PULSE PARAMETERS */
    params.beamRadius        = 0.5e-2;
    params.rabiFrequency     = 2 * pi * 3e6;
    params.psrGradient       = {0., 0., 0.};
    params.laserPhase        = 0.0;

    /* DETECTOR PARAMETERS */
    params.coherenceLength   = 1e-5;

    // create the atom interferometer
    AISAtomInterferometer* ai = new AISAtomInterferometer(params);
    // run, detect, and write
    ai->run();

    // check the status
    writeFile("/home/noammouelle/sim/data_ais++/test.txt", ai->GetAtomEnsemble());

    //ai->detect();
    //ai->write("/home/noammouelle/sim/data_ais++/test.h5");

    // delete the atom interferometer
    //delete ai;

    return 0;
}






