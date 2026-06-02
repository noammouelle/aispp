#ifndef AISATOMENSEMBLE_HH
#define AISATOMENSEMBLE_HH

#include <vector>
#include <memory>

#include "AISConstants.hh"
#include "AISUtilities.hh"
#include "AISAtom.hh"

using atomVector = std::vector<std::unique_ptr<AISAtom>>;

class AISAtomEnsemble
{
private:
    std::unique_ptr<atomVector> fpAtomVector = std::make_unique<atomVector>();

public:
    AISAtomEnsemble(int nAtoms, double transtemperature, double longtemp, double width,
                    doubleThreeVector initialPosition, doubleThreeVector initialVelocity,
                    double seed, int initialState);

    AISAtomEnsemble(double xMin,  double xMax,  int nx,
                    double yMin,  double yMax,  int ny,
                    double zMin,  double zMax,  int nz,
                    double vxMin, double vxMax, int nvx,
                    double vyMin, double vyMax, int nvy,
                    double vzMin, double vzMax, int nvz,
                    int initialState);

    ~AISAtomEnsemble();

    __float128 omega0 = 2 * pi128 * c / 689 * pow(10, 9);

    int GetNumberOfAtoms();

    std::unique_ptr<atomVector>& GetAtomVector();
    std::unique_ptr<AISAtom>& GetAtom(int atomIndex);
};

#endif