#include "AISIo.hh"

void writeAtomEnsembleToFile(std::string fName, AISAtomEnsemble* atomEnsemble)
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