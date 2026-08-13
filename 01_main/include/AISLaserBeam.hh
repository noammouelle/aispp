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

    // Load the sampled beam (datasets x, y, z, phase, amplitude) written by
    // aisoptics' AISPPExporter. Only meaningful for beamType == "interpolated";
    // AISDriver calls it only for that type.
    void SetInterpolationGrids(std::string filename);
    bool HasInterpolationGrids() const;
    // Number of distinct beam files currently held in the shared cache --
    // exposed for tests to pin the sharing behaviour (see InterpolationCache
    // in the private section), not meant for use by simulation code.
    static size_t InterpolationCacheSize();

    void SetTipTilt(double tiptiltX, double tiptiltY);
    std::pair<double, double> GetTipTilt();

private:
    // Shear the sample point to model a tipped/tilted retroreflector.
    doubleThreeVector ApplyTipTilt(const doubleThreeVector& pos) const;
    // Zernike aberration phase over the aperture; 0 when no coefficients are set.
    double GetZernikePhase(const doubleThreeVector& pos) const;

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

    // Sampled beam data (x, y, z, phase, amplitude) loaded by
    // SetInterpolationGrids. Held via a filename-keyed static cache
    // (InterpolationCache()) rather than a per-beam copy: AISDriver builds
    // one AISLaserBeam per pulse, and a real pulse sequence commonly has
    // hundreds to thousands of pulses referencing the same one or two beam
    // files, so a naive per-beam copy multiplies the field's memory
    // footprint by the pulse count -- large enough to exhaust system memory
    // for a big field and a long sequence even though the file itself is
    // only read once.
    struct InterpolationGrids
    {
        std::vector<double> x, y, z, phase, amplitude;
    };
    static std::map<std::string, std::shared_ptr<const InterpolationGrids>>& InterpolationCache();
    std::shared_ptr<const InterpolationGrids> interpolationGrids;

    double tiptiltX = 0.0;
    double tiptiltY = 0.0;

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