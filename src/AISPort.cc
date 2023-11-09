#include "AISPort.hh"

AISPort::AISPort(/* args */)
{
}

AISPort::~AISPort()
{
}

int AISPort::getNumberOfWavePackets()
{
    return wavePacketIndices.size();
}