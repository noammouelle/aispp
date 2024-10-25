#ifndef AISGAUSSIANWAVEFRONT_HH
#define AISGAUSSIANWAVEFRONT_HH

#include "AISWaveFront.hh"
#include "AISUtilities.hh"
#include "AISConstants.hh"

class AISGaussianWaveFront : public AISWaveFront
{
public:
    double beamRadius; // beam radius at the origin of the beam (mirror of lens)
    double w0;
    double zR;

    AISGaussianWaveFront(doubleThreeVector psrGradient, double laserPhase, double beamRadius);
    ~AISGaussianWaveFront();

    double GetNonShiftedValue(const doubleThreeVector& pos, const doubleThreeVector& k);
    double GetValue(const doubleThreeVector& pos, const doubleThreeVector& k) override;

    // relevant functions
    double R(const double& z);
    double w(const double& z);
    double gouyPhase(const double& z);
};

#endif