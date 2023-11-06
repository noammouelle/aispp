#ifndef PORTFRAME_HH
#define PORTFRAME_HH

#include <vector>
#include <algorithm>
#include <cmath>

#include "AISUtilities.hh"
#include "AISPort.hh"
#include "AISAtom.hh"

using portVector = std::vector<AISPort*>;

class AISPortFrame
{
public:
    AISPortFrame(AISAtom* anAtom, intTuple adjacentWavepacketIndices);
    ~AISPortFrame();

    portVector* fpPortVector;

    bool isValueInVector(int value, const intVector& vector);
    int GetNumberOfPorts();
    AISPort* GetPort(int portIndex);
};


#endif