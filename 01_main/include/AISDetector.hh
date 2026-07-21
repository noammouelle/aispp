#ifndef AISDETECTOR_HH
#define AISDETECTOR_HH

#include <array>
#include <vector>
#include <random>

#include "AISUtilities.hh"
#include "AISAtomEnsemble.hh"
#include "AISPortFrame.hh"

using portFrameVector = std::vector<std::unique_ptr<AISPortFrame>>;

class AISDetector
{
public:
    AISDetector(std::unique_ptr<AISAtomEnsemble>& pAtomEnsemble, double coherenceLength, bool ultraFast = false);
    ~AISDetector();

    intVector SampleAllPorts();
    int SamplePort(std::unique_ptr<AISPortFrame>& aPortFrame);

    int GetNumberOfSamples();

//protected:
    intTuple GetAdjacentWavepackets(std::unique_ptr<AISAtom>& anAtom, double coherenceLength);
    doubleVector GetPortProbabilities(int atomIndex, intVector adjacentWavepackets);

    std::unique_ptr<AISPortFrame>& GetPortFrame(int portFrameIndex);
    int GetNumberOfPortFrames();

    double coherenceLength;

    int numSamples = 0;

    std::unique_ptr<portFrameVector> fpPortFrameVector;
};

#endif