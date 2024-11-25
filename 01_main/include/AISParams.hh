#ifndef AISPARAMS_HH
#define AISPARAMS_HH

#include "AISUtilities.hh"
#include <vector>
#include <set>
#include <quadmath.h> // For __float128

class AISParams
{
public:
    AISParams();
    ~AISParams();

    // Cloud parameters
    int nAtoms = 0;
    double cloudRadius = 0.0;
    double cloudTemperature = 0.0;
    doubleThreeVector initialPosition = {0.0, 0.0, 0.0};
    doubleThreeVector initialVelocity = {0.0, 0.0, 0.0};

    // Potential parameters
    std::string potentialType;

    // Laser parameters
    std::vector<double> rabiFrequencies;
    std::vector<double> phi0;
    std::vector<__float128> initialPulseTimes;
    std::vector<__float128> finalPulseTimes;
    std::vector<double> kXVector;
    std::vector<double> kYVector;
    std::vector<double> kZVector;
    std::vector<__float128> omegaVector;
    std::vector<std::string> wavefrontTypeVector;
    std::vector<__float128> betaVector;

    // Kinematic parameters
    __float128 detectionTime;

    // Simulation parameters
    double amplitudeThreshold;
    double coherenceLength;
    bool useMcBranching;

    // Output parameters
    bool printPorts;
};


#endif