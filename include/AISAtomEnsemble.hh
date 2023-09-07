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

    int GetNumberOfAtoms();

    AISAtom* GetAtom(int atomIndex);
};

#endif