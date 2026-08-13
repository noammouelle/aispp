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
    // Angular velocity of the rotating frame, rad/s, in simulation-frame
    // components. Only used by the rotating_* potential types; zero means the
    // frame is inertial. Note the input key is "rotation": "omega" is already
    // taken by the laser angular frequencies.
    doubleThreeVector rotationRate = {0.0, 0.0, 0.0};

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
    std::vector<std::string> beamInterpolationParamsFilenames;
    std::vector<double> tiptiltX;
    std::vector<double> tiptiltY;

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
    bool ultraFast;

    // Output parameters
    bool printPorts;
    bool printWavePackets;
    bool printTrajectory = false;

    // Phase space grid parameters (used when usePhaseSpaceGrid = true)
    bool usePhaseSpaceGrid = false;
    double xGridMin = 0.0,  xGridMax = 0.0;  int nGridX = 1;
    double yGridMin = 0.0,  yGridMax = 0.0;  int nGridY = 1;
    double zGridMin = 0.0,  zGridMax = 0.0;  int nGridZ = 1;
    double vxGridMin = 0.0, vxGridMax = 0.0; int nGridVX = 1;
    double vyGridMin = 0.0, vyGridMax = 0.0; int nGridVY = 1;
    double vzGridMin = 0.0, vzGridMax = 0.0; int nGridVZ = 1;
};


#endif