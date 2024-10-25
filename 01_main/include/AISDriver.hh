#ifndef AISDRIVER_HH
#define AISDRIVER_HH

#include "AISParams.hh"
#include "AISDataIO.hh"
#include "AISAtomEnsemble.hh"
#include "AISDetector.hh"
#include "AISPulsePropagator.hh"
#include "AISKinematicPropagator.hh"
#include "AISLaserBeam.hh"

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
};

#endif