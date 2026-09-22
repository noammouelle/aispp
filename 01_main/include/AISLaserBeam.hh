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

    // Analytic single transverse Fourier mode on the carrier. Inactive when
    // epsilon == 0, so setting it is backward compatible: gaussian and confocal
    // beams behave exactly as before until a mode is supplied. Unlike
    // wtype=interpolated this needs no sampled field and therefore no grid, so
    // there is nothing to fall outside of -- the same call is valid for a
    // source at 0-20 m and one at 100-120 m.
    void SetFourierMode(const FourierMode& mode);
    // wtype "ripple": phi = phi0 + eps cos(q x + psi), grad phi = (-eps q sin(q x + psi), 0, 0), Omega = rabiFreq.
    void SetRipple(double q, double eps, double psi);
    FourierMode GetFourierMode() const;

    // Whether GetDelPhi returns the mode's wavefront gradient or zero. Zero
    // reproduces the historical behaviour (see KNOWN_ISSUES.md), and keeping it
    // switchable is what allows a run with and without the gradient recoil to
    // be compared directly.
    void SetWavefrontGradientEnabled(bool enabled);
    bool GetWavefrontGradientEnabled() const;

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

    FourierMode fourierMode;
    double rippleQ = 0.0, rippleEps = 0.0, ripplePsi = 0.0;   // wtype "ripple" only
    bool wavefrontGradientEnabled = false;  // not yet correct; see KNOWN_ISSUES.md

    // Mode contributions RELATIVE to the carrier, so callers layer them on top
    // of the existing beam-type branches without touching their conventions.
    void FourierModeContributions(const doubleThreeVector& pos,
                                  double& phase, double& amplitudeFactor) const;
    double GetFourierModePhase(const doubleThreeVector& pos) const;
    double GetFourierModeAmplitudeFactor(const doubleThreeVector& pos) const;
    doubleThreeVector GetFourierModeGradPhase(const doubleThreeVector& pos) const;
    // Signed contributions actually applied, honouring beam direction and the
    // per-beam-type phase sign convention.
    double ModePhaseContribution(const doubleThreeVector& pos) const;
    double ModeAmplitudeContribution(const doubleThreeVector& pos) const;
    doubleThreeVector ModeGradContribution(const doubleThreeVector& pos) const;

    // Carrier coordinate for the current beam type, matching GetPhi exactly.
    // Returns the shifted position and the signed distance from the mirror.
    void CarrierCoordinates(const doubleThreeVector& pos,
                            doubleThreeVector& shiftedPos,
                            double& dzFromMirror,
                            double& zSign) const;

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