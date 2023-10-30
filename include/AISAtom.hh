#ifndef ATOM_HH
#define ATOM_HH

#include "AISConstants.hh"
#include <vector>
#include "AISWavePacket.hh"

using wavePacketVector = std::vector<AISWavePacket*>;

class AISAtom{
public:
    AISAtom(doubleThreeVector initialPos, doubleThreeVector initialVel, __float128 initialTime);
    ~AISAtom();

    int GetNumberOfWavePackets();

    AISWavePacket* GetWavePacket(int wavePacketIndex);

    void DeleteWavePackets();
    void AddWavePackets(wavePacketVector* newWavePacketVector);

private:
    //wavePacketVector fWavePacketVector;
    wavePacketVector* fpWavePacketVector = new wavePacketVector;

};

#endif 