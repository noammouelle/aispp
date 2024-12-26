#ifndef AISPORT_HH
#define AISPORT_HH

#include "AISUtilities.hh"

class AISPort
{
public:
    AISPort(/* args */);
    ~AISPort();

    double probabilityAmplitude = 0.0;
    doubleThreeVector position = {0.0, 0.0, 0.0};
    doubleThreeVector velocity = {0.0, 0.0, 0.0};
    double phaseShift = 0.0;
    double phaseShiftError = 0.0;
    int state = 0;
    intVector wavePacketIndices;

    bool interfering;

    int getNumberOfWavePackets();
};


#endif