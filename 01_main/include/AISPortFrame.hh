#ifndef PORTFRAME_HH
#define PORTFRAME_HH

#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>

#include "AISUtilities.hh"
#include "AISPort.hh"
#include "AISAtom.hh"

using portVector = std::vector<std::unique_ptr<AISPort>>;

class AISPortFrame
{
public:
    AISPortFrame(std::unique_ptr<AISAtom>& anAtom, double coherenceLength, bool ultraFast = false);
    ~AISPortFrame();

    double coherenceLength;
    bool ultraFast;

    void initializePortVector(int numberOfPorts);
    void associatePortsWithWavePackets(std::unique_ptr<AISAtom>& anAtom, intTuple adjacentWavepacketIndices);
    void setPortParameters(std::unique_ptr<AISAtom>& anAtom);
    
    std::vector<int> createGroup(std::unique_ptr<AISAtom>& anAtom, int wavePacketIndex, double coherenceLength);

    std::unique_ptr<portVector> fpPortVector;

    int findValueInVector(int value, const intVector& vector);
    bool isValueInVector(int value, const intVector& vector);
    int GetNumberOfPorts();
    std::unique_ptr<AISPort>& GetPort(int portIndex);

    std::vector<std::string> interferingPaths = {};
    std::vector<std::string> GetInterferingPaths();
};


#endif