#include "AISAtom.hh"

AISAtom::AISAtom(doubleThreeVector initialPos, doubleThreeVector initialVel, __float128 initialTime){
    AISWavePacket* initialWavePacket = new AISWavePacket();
    initialWavePacket->SetPosition(initialPos);
    initialWavePacket->SetVelocity(initialVel);
    initialWavePacket->SetTime(initialTime);
    fpWavePacketVector->push_back(initialWavePacket);
};

AISAtom::~AISAtom(){
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
        delete (*fpWavePacketVector)[i];
    }
};

void AISAtom::AddWavePackets(wavePacketVector* newWavePacketVector){
    fpWavePacketVector = newWavePacketVector;
};