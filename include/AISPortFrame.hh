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
    AISPortFrame(AISAtom* anAtom, intTuple adjacentWavepacketIndices);
    ~AISPortFrame();

    void initializePortVector(int numberOfPorts);
    void associatePortsWithWavePackets(AISAtom* anAtom, intTuple adjacentWavepacketIndices);
    void setPortParameters(AISAtom* anAtom, intTuple adjacentWavepacketIndices);

    portVector* fpPortVector;

    int findValueInVector(int value, const intVector& vector);
    bool isValueInVector(int value, const intVector& vector);
    int GetNumberOfPorts();
    AISPort* GetPort(int portIndex);
};


#endif