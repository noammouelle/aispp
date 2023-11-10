#include <array>
#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <chrono>

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

void writePortFile(std::string fName, AISAtomInterferometer* ai)
{
    // write a file to cross check
    std::ofstream outFile(fName);
    std::cout << "Writing to file: " << fName << std::endl;
    outFile << std::fixed; // Use fixed-point notation
    outFile << std::setprecision(15);
    
    outFile << "State, Probability, X, Y, Z, VX, VY, VZ, PhaseShift" << std::endl;

    // get the port frame
    AISDetector* aDetector = ai->fpDetector;
    AISPortFrame* aPortFrame = aDetector->GetPortFrame(0);

    // loop over ports
    for(int portIndex = 0; portIndex < aPortFrame->GetNumberOfPorts(); portIndex++)
    {   
        // get port
        AISPort* port = aPortFrame->GetPort(portIndex);
        // get properties
        int state = port->state;
        double probability = port->probabilityAmplitude;
        doubleThreeVector pos = port->position;
        doubleThreeVector vel = port->velocity;
        double phaseShift = port->phaseShift;
        // write to file
        outFile << state << "," << probability << "," << pos[0] << "," << pos[1] << "," << pos[2]
                << "," << vel[0] << "," << vel[1] << "," << vel[2] << "," << phaseShift << std::endl;

    }

    outFile.close();
}

int main()
{   
    // create a parameter object
    AISAtomInterferometerParams params;

    /* CLOUD PARAMETERS */
    params.nAtoms            = 10000;
    params.cloudTemperature  = 15e-12;
    params.cloudWidth        = 1e-3;

    params.initialPosition   = {0.,0.,0.};
    params.initialVelocity   = {0.,0.,0.};

    /* SEQUENCE PARAMETERS */
    params.interogationTime  = 1.0q;
    params.nSteps            = 20;

    /* PULSE PARAMETERS */
    params.beamRadius        = 1e-2;
    params.rabiFrequency     = 2 * pi * 3e6;
    params.psrGradient       = {1.6e-3 * pi / lambdaSr87, 0., 0.};
    params.laserPhase        = 0.0;

    /* DETECTOR PARAMETERS */
    params.coherenceLength   = 1e-5;

    // create the atom interferometer
    AISAtomInterferometer* ai = new AISAtomInterferometer(params);

    auto t0 = std::chrono::high_resolution_clock::now();
    // run, detect, and write
    ai->run();

    auto t1 = std::chrono::high_resolution_clock::now();
    auto durationRun = std::chrono::duration_cast<std::chrono::seconds>(t1 - t0);
    std::cout<<"Run was successful! ( "<< durationRun.count() << " seconds )" <<std::endl;

    ai->detect();
    auto t2 = std::chrono::high_resolution_clock::now();
    auto durationDetect = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1);
    std::cout << "Detection was successful! ( " << durationDetect.count() << " seconds )" << std::endl;

    ai->write("/home/noammouelle/sim/data_ais++/test.h5");
    auto t3 = std::chrono::high_resolution_clock::now();
    auto durationWrite = std::chrono::duration_cast<std::chrono::seconds>(t3 - t2);
    std::cout << "Write was successful! ( " << durationWrite.count() << " seconds )" << std::endl;

    return 0;
}






