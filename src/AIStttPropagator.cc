#include "AIStttPropagator.hh"

AIStttPropagator::AIStttPropagator(__float128 deltaTime)
{
    fDeltaTime = deltaTime;
    fFreePropagator = new AISFreePropagator();
}

AIStttPropagator::~AIStttPropagator()
{
    delete fFreePropagator;
}

void AIStttPropagator::PropagateEnsemble(AISAtomEnsemble* atomEnsemble)
{
    // propagate half a step freely
    fFreePropagator->PropagateEnsemble(atomEnsemble, fDeltaTime / 2);
    // apply the ttt scheme
    TttPropagateEnsemble(atomEnsemble);
    // propagate half a step freely
    fFreePropagator->PropagateEnsemble(atomEnsemble, fDeltaTime / 2);
}

void AIStttPropagator::TttPropagateEnsemble(AISAtomEnsemble* atomEnsemble)
{
    for(int atomIndex = 0; atomIndex < atomEnsemble->GetNumberOfAtoms() ; ++atomIndex)
    {
        AISAtom* currentAtom = atomEnsemble->GetAtom(atomIndex);
        TttPropagateAtom(currentAtom);
    }
}

void AIStttPropagator::TttPropagateAtom(AISAtom* atom)
{
    for(int wavePacketIndex = 0; wavePacketIndex < atom->GetNumberOfWavePackets(); ++wavePacketIndex)
    {
        // Get pointers to the two output wavepackets
        AISWavePacket* wavePacket1 = atom->GetWavePacket(wavePacketIndex);
        AISWavePacket* wavePacket2 = new AISWavePacket();

        /*
        NEED TO FIND EXPRESSIONS FOR ABCD MATRICES
        */

        // Excited to ground
        for(int dim = 0; dim < 3; ++dim)
        {
            dx[dim] = -GetDoubleWaveVector()[dim]/massSr87 * hbar * fDeltaTime / 2;
            dv[dim] =  GetDoubleWaveVector()[dim]/massSr87 * hbar;
        }
    }
}