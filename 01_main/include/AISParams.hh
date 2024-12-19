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
    int initialState = 0;
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
    std::vector<__float128> frequencyChirpVector;
    std::vector<double> kXVector;
    std::vector<double> kYVector;
    std::vector<double> kZVector;
    std::vector<double> kXChirpVector, kYChirpVector, kZChirpVector;
    std::vector<__float128> omegaVector;
    std::vector<std::string> wavefrontTypeVector;

    // Kinematic parameters
    __float128 detectionTime;

    // Simulation parameters
    double amplitudeThreshold;
    double coherenceLength;
    double seed;
    bool useMcBranching;
    bool ignoreDetuning;
    double xDetMin, xDetMax, yDetMin, yDetMax, zDetMin, zDetMax;
    bool useDetVolSelection;

    // Output parameters
    bool printPorts;
    bool printWavePackets;
};


#endif