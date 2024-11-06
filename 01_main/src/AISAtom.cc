#include "AISAtom.hh"

AISAtom::AISAtom(doubleThreeVector initialPos, doubleThreeVector initialVel, __float128 initialTime)
{
    std::unique_ptr<AISWavePacket> initialWavePacket(new AISWavePacket());
    initialWavePacket->SetPosition(initialPos);
    initialWavePacket->SetVelocity(initialVel);
    initialWavePacket->SetTime(initialTime);

    // init wavepacket vector
    fpWavePacketVector = std::make_unique<wavePacketVector>();
    // add the initial wavepacket to the vector
    fpWavePacketVector->push_back(std::move(initialWavePacket));
};

AISAtom::~AISAtom(){
    // free all the wavepacket smart pointers and then free the vector smart pointer
    fpWavePacketVector->clear();
    fpWavePacketVector.reset();
};

int AISAtom::GetNumberOfWavePackets(){
    return fpWavePacketVector->size();
};

std::unique_ptr<AISWavePacket>& AISAtom::GetWavePacket(int wavePacketIndex)
{
    return fpWavePacketVector->at(wavePacketIndex);
}

void AISAtom::DeleteWavePackets(){
    fpWavePacketVector->clear();
    fpWavePacketVector.reset();
};

void AISAtom::AddWavePackets(std::unique_ptr<wavePacketVector>& newWavePacketVector){
    fpWavePacketVector = std::move(newWavePacketVector);
};