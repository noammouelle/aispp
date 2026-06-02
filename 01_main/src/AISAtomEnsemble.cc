#include <cmath>
#include <random>

#include "AISAtomEnsemble.hh"

AISAtomEnsemble::AISAtomEnsemble(int nAtoms, double transtemperature, double longtemperature, double width,
                                 doubleThreeVector initialPosition, doubleThreeVector initialVelocity,
                                 double seed, int initialState)
{
    // Isotropic Gaussian distrib, can be changed later by
    // implementing a generator class
    std::random_device rd;
    std::mt19937 gen;
    if (seed == -1) {
        gen.seed(rd());
    } else {
        gen.seed(seed);
    }
    doubleThreeVector centralPos = initialPosition;
    double stdPos = width;
    doubleThreeVector centralVel = initialVelocity;
    double stdVelZ = sqrt(longtemperature * kB / massSr87);
    double stdVelX = sqrt(transtemperature * kB / massSr87);

    std::normal_distribution<double> posDistributionX(centralPos[0], stdPos); // X
    std::normal_distribution<double> velDistributionX(centralVel[0], stdVelX);

    std::normal_distribution<double> posDistributionY(centralPos[1], stdPos); // Y
    std::normal_distribution<double> velDistributionY(centralVel[1], stdVelX);
    
    std::normal_distribution<double> posDistributionZ(centralPos[2], stdPos); // Z
    std::normal_distribution<double> velDistributionZ(centralVel[2], stdVelZ);

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
        std::unique_ptr<AISAtom> currentAtom(new AISAtom(sampledPos[atom_i],
                                                         sampledVel[atom_i],
                                                         0.0q, initialState));
        fpAtomVector->push_back(std::move(currentAtom));
    }
};

AISAtomEnsemble::AISAtomEnsemble(double xMin,  double xMax,  int nx,
                                 double yMin,  double yMax,  int ny,
                                 double zMin,  double zMax,  int nz,
                                 double vxMin, double vxMax, int nvx,
                                 double vyMin, double vyMax, int nvy,
                                 double vzMin, double vzMax, int nvz,
                                 int initialState)
{
    auto linspace = [](double lo, double hi, int n) {
        std::vector<double> v(n);
        if (n == 1) { v[0] = lo; return v; }
        for (int i = 0; i < n; ++i)
            v[i] = lo + i * (hi - lo) / (n - 1);
        return v;
    };

    auto xs  = linspace(xMin,  xMax,  nx);
    auto ys  = linspace(yMin,  yMax,  ny);
    auto zs  = linspace(zMin,  zMax,  nz);
    auto vxs = linspace(vxMin, vxMax, nvx);
    auto vys = linspace(vyMin, vyMax, nvy);
    auto vzs = linspace(vzMin, vzMax, nvz);

    for (double x : xs)
    for (double y : ys)
    for (double z : zs)
    for (double vx : vxs)
    for (double vy : vys)
    for (double vz : vzs)
    {
        doubleThreeVector pos = {x, y, z};
        doubleThreeVector vel = {vx, vy, vz};
        fpAtomVector->push_back(std::make_unique<AISAtom>(pos, vel, 0.0q, initialState));
    }
};

AISAtomEnsemble::~AISAtomEnsemble()
{
    // free all the unique_ptr to the atoms then free the unique_ptr to the atom vector
    fpAtomVector->clear();
    fpAtomVector.reset();
};

int AISAtomEnsemble::GetNumberOfAtoms()
{
    return fpAtomVector->size();
}

std::unique_ptr<atomVector>& AISAtomEnsemble::GetAtomVector()
{
    return fpAtomVector;
}

std::unique_ptr<AISAtom>& AISAtomEnsemble::GetAtom(int atomIndex)
{
    return fpAtomVector->at(atomIndex); // Find better way
}
