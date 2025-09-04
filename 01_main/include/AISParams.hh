#ifndef AISPARAMS_HH
#define AISPARAMS_HH

#include "AISUtilities.hh"
#include <vector>
#include <set>
#include <quadmath.h> // For __float128
#include <map>

class AISParams
{
public:
    AISParams();
    ~AISParams();

    // Cloud parameters
    int initialState = 0;
    int nAtoms = 0;
    double cloudRadius = 0.0;
    double cloudTransTemperature = 0.0;
    double cloudLongTemperature = 0.0;
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
    std::vector<double> waistVector;
    std::vector<double> zLaserVector; // z position of the laser source
    std::vector<double> focalLengthVector; // focal length of the lenses
    std::map<int, std::vector<double>> zernikeCoeff; // index -> coeff map
    std::vector<double> beamRadiusVector;
    std::vector<double> baselineLengthVector;

    // Kinematic parameters
    __float128 detectionTime;

    // Simulation parameters
    double amplitudeThreshold;
    double coherenceLength;
    double seed;
    bool useMcBranching;
    bool ignoreDetuning;
    bool usePathSelection;
    bool useStaticApprox;
    std::vector<std::string> pathsToSimulate;
    double xDetMin, xDetMax, yDetMin, yDetMax, zDetMin, zDetMax;
    bool useDetVolSelection;
    double gslQagAbsError, gslQagRelError, gslKinAbsError, gslKinRelError, gslPulseAbsError, gslPulseRelError;

    // Output parameters
    bool printPorts;
    bool printWavePackets;
};


#endif