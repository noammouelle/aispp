#ifndef AISLASERBEAM_HH
#define AISLASERBEAM_HH

#include "AISConstants.hh"
#include "AISUtilities.hh"
#include "AISWavefronts.hh"
#include "AISEnvelopes.hh"

#include <memory>
#include <map>

#include <gsl/gsl_deriv.h>
#include <H5Cpp.h>

class AISLaserBeam
{
public:
    AISLaserBeam(doubleThreeVector k, __float128 omega, double rabiFreq, double phi0);
    ~AISLaserBeam();

    // At the moment I suppose that I can somehow pass a function as arguments
    // of the constructor to get the phases etc. Need to check on WiFi
    double GetPhi(const doubleThreeVector& pos); 
    doubleThreeVector GetDelPhi(const doubleThreeVector& pos);
    double GetRabiFreq(const doubleThreeVector& pos, const __float128& t0, const __float128& t);
    
    doubleThreeVector GetK(__float128 t);
    void SetK(doubleThreeVector k);

    __float128 GetOmega(__float128 t);
    void SetOmega(__float128 omega);

    __float128 GetFrequencyChirp();
    void SetFrequencyChirp(__float128 frequencyChirp);

    doubleThreeVector GetKChirp();
    void SetKChirp(doubleThreeVector kChirp);

    double GetW0();
    void SetW0(double w0);

    double GetZLaser();
    void SetZLaser(double zLaser);

    double GetFocalLength();
    void SetFocalLength(double focalLength);

    void SetZernikeCoeffs(std::map<int, double> zernikeCoeffs);
    std::map<int, double> GetZernikeCoeffs();

    void SetBeamRadius(double beamRadiusValue);
    double GetBeamRadius();

    void SetBeamType(std::string beamTypeValue);
    std::string GetBeamType();

    void SetBaselineLength(double baselineLengthValue);
    double GetBaselineLength();

    static double GetPhiWrapper(double xi, void* params);

    void SetInterpolationGrids(std::string filename);

    void SetTipTilt(double tiptiltX, double tiptiltY);
    std::pair<double, double> GetTipTilt();

private:
    double rabiFreq; // central rabi frequency

    doubleThreeVector k;
    __float128 omega, frequencyChirp;
    doubleThreeVector kChirp;
    double w0;
    double phi0;
    double zLaser;
    double focalLength;
    double beamRadius;
    double baselineLength;

    std::map<int, double> zernikeCoeffs;

    std::vector<double> interpolationGridX;
    std::vector<double> interpolationGridY;
    std::vector<double> interpolationGridZ;
    std::vector<double> phaseInterpolationGridValues;
    std::vector<double> amplitudeInterpolationGridValues;

    double tiptiltX;
    double tiptiltY;

    // options params
    std::string beamType;

    //std::shared_ptr<wavefrontFunctionType> wavefrontFunction;
    //std::shared_ptr<delWavefrontFunctionType> delWavefrontFunction;
    //std::shared_ptr<rabifreqFunctionType> rabifreqFunction;
};

struct getPhiWrapperParams
{
    doubleThreeVector pos;
    int index;
    std::shared_ptr<AISLaserBeam> beam;
};


#endif