#include "AISPulsePropagator.hh"

AISPulsePropagator::AISPulsePropagator(AISLaserBeam* beam, __float128 t0, __float128 t1)
{
    this->laserBeam=beam;
    this->t0=t0;
    this->t1=t1;
}

AISPulsePropagator::~AISPulsePropagator()
{
    delete this->laserBeam;
}

void AISPulsePropagator::PropagateEnsemble(AISAtomEnsemble* atomEnsemble)
{
    #pragma omp parallel for
    for(int i_atom = 0; i_atom < atomEnsemble->GetNumberOfAtoms(); ++i_atom)
    {
        AISAtom* currentAtom = atomEnsemble->GetAtom(i_atom);
        PropagateAtom(currentAtom);
    }
}

void AISPulsePropagator::PropagateAtom(AISAtom* atom)
{
    for(int i_atom)
}