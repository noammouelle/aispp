#ifndef AISATOMINTERFEROMETER_HH
#define AISATOMINTERFEROMETER_HH

#include <string>
#include <H5Cpp.h>

#include "AISAtomEnsemble.hh"
#include "AISAtomInterferometerParams.hh"
#include "AISLinearGravityPropagator.hh"
#include "AIStttPropagator.hh"
#include "AISDetector.hh"
#include "AISUtilities.hh"

class AISAtomInterferometer
{
public:
    AISAtomInterferometer(AISAtomInterferometerParams params);
    ~AISAtomInterferometer();

    void run();
    void detect();
    void write(std::string filename);

    AISAtomEnsemble* GetAtomEnsemble();

//private:
    AISAtomEnsemble* atomEnsemble;
    AIStttPropagator* beamSplitterPropagator1;
    AIStttPropagator* beamSplitterPropagator2;
    AIStttPropagator* mirrorPropagator;
    AISLinearGravityPropagator* driftPropagator;
    AISLinearGravityPropagator* interrogationTimePropagator;

    AISAtomInterferometerParams fParams;

    AISDetector* fpDetector;
    double coherenceLength;

    doubleThreeVector computeDetunedWaveVector(const doubleThreeVector& k0, const double& vz);
    __float128 computeDetunedOmega(const __float128& omega0, const double& vz);

    intVector sampledPortIndices;
};



#endif