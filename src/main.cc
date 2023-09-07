#include <array>
#include <vector>
#include <fstream>

#include "AISAtomEnsemble.hh"
#include "AISFreePropagator.hh"
#include "AISLinearGravityPropagator.hh"
#include "AISUtilities.hh"

void writeFile(std::string fName, AISAtomEnsemble* atomEnsemble)
{
    // write a file to cross check
    std::ofstream outFile(fName);
    outFile << "X, Y, Z, VX, VY, VZ" << std::endl;
    for(int i = 0; i < atomEnsemble->GetNumberOfAtoms(); ++i)
    {
        AISAtom* currentAtom   = atomEnsemble->GetAtom(i);
        doubleThreeVector currentPos = currentAtom->GetWavePacket(0)->GetPosition();
        doubleThreeVector currentVel = currentAtom->GetWavePacket(0)->GetVelocity();

        outFile << currentPos[0] << "," << currentPos[1] << "," << currentPos[2] << "," 
                << currentVel[0] << "," << currentVel[1] << "," << currentVel[2] << std::endl;
    }

    outFile.close();
}

int main()
{
    int numAtoms      = 100000;
    double cloudTemp  = 15e-9;
    double cloudWidth = 1e-3;

    doubleThreeVector initialPos = {0.,0.,0.};
    doubleThreeVector initialVel = {0.,0.,0.};

    AISAtomEnsemble* atomEnsemble = new AISAtomEnsemble(numAtoms, cloudTemp, cloudWidth,
                                                        initialPos, initialVel);
    
    AISLinearGravityPropagator* propagator = new AISLinearGravityPropagator();

    writeFile("/home/noammouelle/sim/ais++/output0.csv",atomEnsemble);

    propagator->PropagateEnsemble(atomEnsemble,1.0q);

    writeFile("/home/noammouelle/sim/ais++/output1.csv",atomEnsemble);

    return 0;
}






