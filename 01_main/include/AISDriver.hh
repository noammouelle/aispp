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
// time
#include <chrono>

class AISDriver
{
public:
    AISDriver(AISParams params);
    ~AISDriver();

    void Run();
    void RunPathFinder();
    void RunAll();
    void Detect();
    void WriteDetectedAtomsToFile(std::string fName);
    void WriteWavePacketsToFile(std::string fName);
    void WritePortsToFile(std::string fName);
    void WritePhaseSpaceMapToFile(std::string fName);
    void WriteTrajectoryToFile(std::string fName);

    void RecordSnapshot(double time, const std::string& label);

    // Write the frame angular velocity (rad/s) as a "rotation" dataset, so that
    // downstream analysis can tell a rotating run from an inertial one.
    void WriteRotationRate(H5::H5File& file);

    std::vector<std::string> GetDetectablePaths(doubleThreeVector pos0, doubleThreeVector vel0);

// private:
    AISParams params;
    std::unique_ptr<AISAtomEnsemble> atomEnsemble;
    std::unique_ptr<AISDetector> detector;
    std::shared_ptr<AISKinematicPropagator> kinematicPropagator;

    std::vector<std::shared_ptr<AISPulsePropagator>> pulsePropagators;
    std::vector<std::unique_ptr<AISLaserBeam>> laserBeams;

    std::vector<int> sampledPortIndices;
    std::vector<std::string> interferingPaths;
    std::map<std::string, doubleThreeVector> pathToFinalPosMap;

    // Trajectory snapshots (populated when printTrajectory = true)
    std::vector<double>      fTrajTimes;
    std::vector<std::string> fTrajLabels;
    std::vector<int>         fTrajSnapIdx;
    std::vector<int>         fTrajAtomIdx;
    std::vector<std::string> fTrajPaths;
    std::vector<int>         fTrajStates;
    std::vector<double>      fTrajAmplitudes;
    std::vector<std::array<double,3>> fTrajPositions;
    std::vector<std::array<double,3>> fTrajVelocities;
};

#endif