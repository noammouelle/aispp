#ifndef PORTFRAME_HH
#define PORTFRAME_HH

#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>

#include "AISUtilities.hh"
#include "AISPort.hh"
#include "AISAtom.hh"

using portVector = std::vector<AISPort*>;

class AISPortFrame
{
public:
    AISPortFrame(AISAtom* anAtom, double coherenceLength);
    ~AISPortFrame();

    double coherenceLength;

    void initializePortVector(int numberOfPorts);
    void associatePortsWithWavePackets(AISAtom* anAtom, intTuple adjacentWavepacketIndices);
    void setPortParameters(AISAtom* anAtom);
    
    std::vector<int> createGroup(AISAtom* anAtom, int wavePacketIndex, double coherenceLength);

    portVector* fpPortVector;

    int findValueInVector(int value, const intVector& vector);
    bool isValueInVector(int value, const intVector& vector);
    int GetNumberOfPorts();
    AISPort* GetPort(int portIndex);
};


#endif