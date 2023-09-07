#ifndef ATOM_HH
#define ATOM_HH

#include <vector>
#include "AISWavePacket.hh"

using wavePacketVector = std::vector<AISWavePacket*>;

class AISAtom{
public:
    AISAtom(doubleThreeVector initialPos, doubleThreeVector initialVel, __float128 initialTime);
    ~AISAtom();

    int GetNumberOfWavePackets();

    AISWavePacket* GetWavePacket(int wavePacketIndex);

private:
    wavePacketVector fWavePacketVector;
    wavePacketVector* fpWavePacketVector = &fWavePacketVector;

};

#endif 