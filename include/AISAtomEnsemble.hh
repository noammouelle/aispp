#ifndef AISATOMENSEMBLE_HH
#define AISATOMENSEMBLE_HH

#include <vector>

#include "AISAtom.hh"

using atomVector = std::vector<AISAtom*>;

class AISAtomEnsemble
{
private:
    atomVector fAtomVector;
    atomVector* fpAtomVector = &fAtomVector;
public:
    AISAtomEnsemble(int nAtoms, double temperature, double width,
                    doubleThreeVector initialPosition, doubleThreeVector initialVelocity);
    ~AISAtomEnsemble();

    __float128 omega0 = 2 * pi128 * c / 689 * pow(10, 9);

    int GetNumberOfAtoms();

    AISAtom* GetAtom(int atomIndex);
};

#endif