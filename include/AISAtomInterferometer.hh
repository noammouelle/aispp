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
#include "AISLmtTttPropagator.hh"

class AISAtomInterferometer
{
public:
    AISAtomInterferometer(AISAtomInterferometerParams params);
    ~AISAtomInterferometer();

    void run();
    void detect();
    void writeDetectedAtomsInfo(std::string filename);
    void writeWavepacketsInfo(std::string filename);

    AISAtomEnsemble* GetAtomEnsemble();

//private:
    AISAtomEnsemble* atomEnsemble;

    AIStttPropagator* beamSplitterPropagator1;
    AIStttPropagator* beamSplitterPropagator2;
    AIStttPropagator* mirrorPropagator;

    AISLmtTttPropagator* lmtPropagator1;
    AISLmtTttPropagator* lmtPropagator2;
    AISLmtTttPropagator* lmtPropagator3;
    AISLmtTttPropagator* lmtPropagator4;
    
    AISLinearGravityPropagator* driftPropagator;
    AISLinearGravityPropagator* interrogationTimePropagator;

    AISAtomInterferometerParams fParams;

    AISDetector* fpDetector;
    double coherenceLength;

    doubleThreeVector computeDetunedWaveVector(const doubleThreeVector& k0, const double& vz);
    __float128 computeDetunedOmega(const __float128& omega0, const double& vz);

    void runAtom(AISAtom* atom);

    intVector sampledPortIndices;
};



#endif