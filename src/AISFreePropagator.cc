#include "AISFreePropagator.hh"

AISFreePropagator::AISFreePropagator(__float128 dt)
{
    deltaTime = dt;
}

AISFreePropagator::~AISFreePropagator()
{}

void AISFreePropagator::SetAddEnergyPhase(bool addEnergyPhase)
{
    fAddEnergyPhase = addEnergyPhase;
}

void AISFreePropagator::PropagateEnsemble(AISAtomEnsemble* atomEnsemble)
{
    #pragma omp parallel for
    for(int i_atom = 0; i_atom < atomEnsemble->GetNumberOfAtoms(); ++i_atom)
    {
        AISAtom* currentAtom = atomEnsemble->GetAtom(i_atom);
        PropagateAtom(currentAtom);
    }
}

void AISFreePropagator::PropagateAtom(AISAtom* atom)
{
    for(int i_wavePacket = 0; i_wavePacket < atom->GetNumberOfWavePackets(); ++i_wavePacket)
    {
        AISWavePacket* currentWavePacket = atom->GetWavePacket(i_wavePacket);
        PropagateWavePacket(currentWavePacket);
    }
}

void AISFreePropagator::PropagateWavePacket(AISWavePacket* wavePacket)
{
    doubleThreeVector currentPos  = wavePacket->GetPosition();
    doubleThreeVector currentVel  = wavePacket->GetVelocity();

    double currentPhaseDouble = wavePacket->GetPhaseDouble();

    __float128 currentTime = wavePacket->GetTime();

    wavePacket->SetPosition(CalculateNewPos(currentPos, currentVel));
    wavePacket->SetVelocity(CalculateNewVel(currentPos, currentVel));
    wavePacket->SetPhaseDouble(CalculateNewPhaseDouble(currentPhaseDouble, 
                                                        currentPos, wavePacket->GetPosition(), 
                                                        currentVel, wavePacket->GetVelocity()));
    if(fAddEnergyPhase && wavePacket->GetState() == 1)
    {
        wavePacket->SetPhaseQuad(CalculateNewPhaseQuad(wavePacket->GetPhaseQuad()));
    };
    
    wavePacket->SetTime(currentTime + deltaTime);
}
