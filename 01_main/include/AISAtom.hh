#ifndef ATOM_HH
#define ATOM_HH

#include "AISConstants.hh"
#include <vector>
#include "AISWavePacket.hh"
#include "AISUtilities.hh"

using wavePacketVector = std::vector<std::unique_ptr<AISWavePacket>>;

class AISAtom{
public:
    AISAtom(doubleThreeVector initialPos, doubleThreeVector initialVel, __float128 initialTime);
    ~AISAtom();

    int GetNumberOfWavePackets();

    AISWavePacket* GetWavePacket(int wavePacketIndex);

    void DeleteWavePackets();
    void AddWavePackets(std::unique_ptr<wavePacketVector> newWavePacketVector);

private:
    std::unique_ptr<wavePacketVector> fpWavePacketVector;

};

#endif 