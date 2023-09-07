#include "AISFreePropagator.hh"

AISFreePropagator::AISFreePropagator()
{}

AISFreePropagator::~AISFreePropagator()
{}

void AISFreePropagator::PropagateEnsemble(AISAtomEnsemble* atomEnsemble, __float128 time)
{
    for(int i_atom = 0; i_atom < atomEnsemble->GetNumberOfAtoms(); ++i_atom)
    {
        AISAtom* currentAtom = atomEnsemble->GetAtom(i_atom);
        PropagateAtom(currentAtom, time);
    }
}

void AISFreePropagator::PropagateAtom(AISAtom* atom, __float128 time)
{
    for(int i_wavePacket = 0; i_wavePacket < atom->GetNumberOfWavePackets(); ++i_wavePacket)
    {
        AISWavePacket* currentWavePacket = atom->GetWavePacket(i_wavePacket);
        PropagateWavePacket(currentWavePacket, time);
    }
}

void AISFreePropagator::PropagateWavePacket(AISWavePacket* wavePacket, __float128 time)
{
    currentPos  = wavePacket->GetPosition();
    currentVel  = wavePacket->GetVelocity();
    currentTime = wavePacket->GetTime();

    wavePacket->SetPosition(CalculateNewPos(currentPos, currentVel, time));
    wavePacket->SetVelocity(CalculateNewVel(currentPos, currentVel, time));
}
