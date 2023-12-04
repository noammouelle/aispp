#include <cmath>
#include <random>

#include "AISConstants.hh"
#include "AISAtomEnsemble.hh"

AISAtomEnsemble::AISAtomEnsemble(int nAtoms, double temperature, double width,
                                 doubleThreeVector initialPosition, doubleThreeVector initialVelocity)
{
    // Isotropic Gaussian distrib, can be changed later by
    // implementing a generator class

    doubleThreeVector centralPos = initialPosition;
    double stdPos = width;
    doubleThreeVector centralVel = initialVelocity;
    double stdVel = sqrt(temperature * kB / massSr87);

    std::random_device rd;
    std::mt19937 gen(rd());

    std::normal_distribution<double> posDistributionX(centralPos[0], stdPos); // X
    std::normal_distribution<double> velDistributionX(centralVel[0], stdVel);

    std::normal_distribution<double> posDistributionY(centralPos[1], stdPos); // Y
    std::normal_distribution<double> velDistributionY(centralVel[1], stdVel);
    
    std::normal_distribution<double> posDistributionZ(centralPos[2], stdPos); // Z
    std::normal_distribution<double> velDistributionZ(centralVel[2], stdVel);

    std::vector<doubleThreeVector> sampledPos(nAtoms);
    std::vector<doubleThreeVector> sampledVel(nAtoms);

    for(int i_sample = 0; i_sample < nAtoms; ++i_sample){

        sampledPos[i_sample][0] = posDistributionX(gen);
        sampledVel[i_sample][0] = velDistributionX(gen);

        sampledPos[i_sample][1] = posDistributionY(gen);
        sampledVel[i_sample][1] = velDistributionY(gen);

        sampledPos[i_sample][2] = posDistributionZ(gen);
        sampledVel[i_sample][2] = velDistributionZ(gen);

    }

    for(int atom_i = 0; atom_i < nAtoms; ++atom_i){
        AISAtom* currentAtom = new AISAtom(sampledPos[atom_i],
                                           sampledVel[atom_i],
                                           0.0q);
        
        fpAtomVector->push_back(currentAtom);
    }
};

AISAtomEnsemble::~AISAtomEnsemble()
{
    delete fpAtomVector;
};

int AISAtomEnsemble::GetNumberOfAtoms()
{
    return fpAtomVector->size();
}

atomVector* AISAtomEnsemble::GetAtomVector()
{
    return fpAtomVector;
}

AISAtom* AISAtomEnsemble::GetAtom(int atomIndex)
{
    return fpAtomVector->at(atomIndex); // Find better way
}
