#include <array>
#include <vector>
#include <fstream>
#include <iomanip>

#include "AISAtomEnsemble.hh"
#include "AISFreePropagator.hh"
#include "AISLinearGravityPropagator.hh"
#include "AISUtilities.hh"

void writeFile(std::string fName, AISAtomEnsemble* atomEnsemble)
{
    // write a file to cross check
    std::ofstream outFile(fName);

    outFile << std::fixed; // Use fixed-point notation
    outFile << std::setprecision(15);

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
    int numAtoms      = 1000;
    double cloudTemp  = 15e-9;
    double cloudWidth = 1e-3;

    doubleThreeVector initialPos = {0.,0.,0.};
    doubleThreeVector initialVel = {0.,0.,0.};

    AISAtomEnsemble* atomEnsemble = new AISAtomEnsemble(numAtoms, cloudTemp, cloudWidth,
                                                        initialPos, initialVel);
    
    __float128 dt = 0.1q;
    int nsteps = 10;

    AISLinearGravityPropagator* propagator = new AISLinearGravityPropagator(dt);

    writeFile("/home/noammouelle/sim/data_ais++/output0.csv",atomEnsemble);

    for(int i = 0; i < nsteps; ++i){
        propagator->PropagateEnsemble(atomEnsemble,dt);
        if(i%10 == 0){
            std::cout << i << "/" << nsteps << std::endl;
        }
    }

    writeFile("/home/noammouelle/sim/data_ais++/output1.csv",atomEnsemble);
    
    return 0;
}






