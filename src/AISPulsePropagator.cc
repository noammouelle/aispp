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
    // create a new wavepacket vector
    wavePacketVector* newWavePackets = new wavePacketVector;

    // Step 1: Û3 * Û2 * Û1 (forward transformations)
    for(int wavePacketIndex = 0; wavePacketIndex < atom->GetNumberOfWavePackets(); ++wavePacketIndex)
    {
        // Get the pointer to the wavepackets
        AISWavePacket* wavepacket0 = atom->GetWavePacket(wavePacketIndex);
        AISWavePacket* wavepacket1 = new AISWavePacket();

        // Apply U1(t_0,t_0) transformation (do nothing)
        // Apply U2(t_0,t_0) transformation
        ApplyU2(wavepacket0, t0, t0);
        // Apply U3(t_0,t)
        ApplyU3(wavepacket0, wavepacket1, t0, t1);
        
        // add the wavepackets to the vector
        newWavePackets->push_back(wavepacket0);
        newWavePackets->push_back(wavepacket1);        
    }

    // Step 2: update the wavepacket vector
    atom->DeleteWavePackets();
    atom->AddWavePackets(newWavePackets);

    // Step 3: perform the inverse transformations
    // Û1^dagger * Û_2^dagger
    for(int wavePacketIndex = 0; wavePacketIndex < atom->GetNumberOfWavePackets(); ++wavePacketIndex)
    {
        AISWavePacket* currentWavepacket = atom->GetWavePacket(wavePacketIndex);
        ApplyU2Dagger(currentWavepacket, t0, t1);
        ApplyU1(currentWavepacket, t1, t0);
    }
}

void AISPulsePropagator::ApplyU1(AISWavePacket* wavepacket, __float128 t0, __float128 t1)
{
    // U1 is just kinematic propagation from t1 to t0 (note the reverse time order)
    
}