#ifndef AISDRIVER_HH
#define AISDRIVER_HH

#include "AISParams.hh"
#include "AISDataIO.hh"
#include "AISAtomEnsemble.hh"
#include "AISDetector.hh"
#include "AISPulsePropagator.hh"
#include "AISKinematicPropagator.hh"
#include "AISLaserBeam.hh"

#include "AISPotentials.hh"
#include "AISWavefronts.hh"
#include "AISEnvelopes.hh"

//hdf5
#include "H5Cpp.h"

class AISDriver
{
public:
    AISDriver(AISParams params);
    ~AISDriver();

    void Run();
    void Detect();
    void WriteDetectedAtomsToFile(std::string fName);
    void WriteWavePacketsToFile(std::string fName);

private:
    AISParams params;
    AISAtomEnsemble* atomEnsemble;
    AISDetector* detector;
    AISKinematicPropagator* kinematicPropagator;

    std::vector<AISPulsePropagator*> pulsePropagators;
    std::vector<AISLaserBeam*> laserBeams;

    std::vector<int> sampledPortIndices;
};

#endif