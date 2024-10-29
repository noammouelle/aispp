#include "AISAtom.hh"

AISAtom::AISAtom(doubleThreeVector initialPos, doubleThreeVector initialVel, __float128 initialTime){
    AISWavePacket* initialWavePacket = new AISWavePacket();
    initialWavePacket->SetPosition(initialPos);
    initialWavePacket->SetVelocity(initialVel);
    initialWavePacket->SetTime(initialTime);
    fpWavePacketVector->push_back(initialWavePacket);
};

AISAtom::~AISAtom(){
    // delete all the wavepackets
    for(int wavePacketIndex = 0; wavePacketIndex < fpWavePacketVector->size(); wavePacketIndex++)
    {
        //delete (*fpWavePacketVector)[wavePacketIndex];
        delete fpWavePacketVector->at(wavePacketIndex);
    }
    // delete the vector
    delete fpWavePacketVector;
};

int AISAtom::GetNumberOfWavePackets(){
    return fpWavePacketVector->size();
};

AISWavePacket* AISAtom::GetWavePacket(int wavePacketIndex)
{
    return (*fpWavePacketVector)[wavePacketIndex];
}

void AISAtom::DeleteWavePackets(){
    for(int i = 0; i < fpWavePacketVector->size(); ++i){
        std::cout << "Deleting wavepacket " << i << std::endl;
        delete (*fpWavePacketVector)[i];
    }
};

void AISAtom::AddWavePackets(wavePacketVector* newWavePacketVector){
    // if the length of the vector is non-zero, add the wavepackets to the atom
    if(newWavePacketVector->size() > 0){
        fpWavePacketVector = newWavePacketVector;
    }
};