#ifndef AISDETECTOR_HH
#define AISDETECTOR_HH

#include <array>
#include <vector>
#include <random>

#include "AISUtilities.hh"
#include "AISAtomEnsemble.hh"
#include "AISPortFrame.hh"

using portFrameVector = std::vector<AISPortFrame*>;

class AISDetector
{
public:
    AISDetector(AISAtomEnsemble* pAtomEnsemble, double coherenceLength);
    ~AISDetector();

    intVector SampleAllPorts();
    int SamplePort(AISPortFrame* aPortFrame);

    int GetNumberOfSamples();

//protected:
    intTuple GetAdjacentWavepackets(AISAtom* anAtom, double coherenceLength);
    doubleVector GetPortProbabilities(int atomIndex, intVector adjacentWavepackets);

    AISPortFrame* GetPortFrame(int portFrameIndex);
    int GetNumberOfPortFrames();

    double coherenceLength;

    int numSamples = 0;

    portFrameVector* fpPortFrameVector;
};

#endif