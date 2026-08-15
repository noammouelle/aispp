#include "AISLaserBeam.hh"

AISLaserBeam::AISLaserBeam(doubleThreeVector aK, __float128 aOmega, double aRabiFreq, double phi0) : k(aK), omega(aOmega), rabiFreq(aRabiFreq), phi0(phi0)
{}

AISLaserBeam::~AISLaserBeam()
{}

doubleThreeVector AISLaserBeam::GetK(__float128 t)
{
    return matrixAdd(this->k,scalarMultiply(this->kChirp,0.5 * static_cast<double>(t)));
}

void AISLaserBeam::SetK(doubleThreeVector kValue)
{
    k = kValue;
}

doubleThreeVector AISLaserBeam::GetKChirp()
{
    return kChirp;
}
void AISLaserBeam::SetKChirp(doubleThreeVector kChirpValue)
{
    kChirp = kChirpValue;
}

/**
 * @brief Calculates the laser frequency at a given time.
 *
 * This function computes the angular frequency (omega) at a specific time (t)
 * by adding the initial angular frequency (omega) to the product of the 
 * frequency chirp rate and the time.
 *
 * @param t The time at which to calculate the angular frequency.
 * @return The angular frequency at time t.
 */
__float128 AISLaserBeam::GetOmega(__float128 t)
{
    return omega + 0.5q * frequencyChirp * t;
}

void AISLaserBeam::SetOmega(__float128 omegaValue)
{
    this->omega=omegaValue;
}

__float128 AISLaserBeam::GetFrequencyChirp()
{
    return frequencyChirp;
}
void AISLaserBeam::SetFrequencyChirp(__float128 frequencyChirpValue)
{
    frequencyChirp=frequencyChirpValue;
}

doubleThreeVector AISLaserBeam::ApplyTipTilt(const doubleThreeVector& pos) const
{
    // A retroreflector tipped by theta shears the beam frame by tan(theta) per
    // unit z. Identity when both angles are zero, so untilted runs are unchanged.
    if(tiptiltX == 0.0 && tiptiltY == 0.0)
        return pos;
    return {pos[0] - pos[2] * tan(tiptiltX * (pi/180.0)),
            pos[1] - pos[2] * tan(tiptiltY * (pi/180.0)),
            pos[2]};
}

double AISLaserBeam::GetZernikePhase(const doubleThreeVector& pos) const
{
    // Wavefront aberration expanded in Zernike polynomials over the aperture,
    // phi = sum_N c_N Z_N(r/beamRadius). Restored in v0.0.2: this block was
    // dropped from GetPhi by commit adbebdf ("Working version"), which silently
    // made every zernikecoeff_N input inert. The loop is empty when no
    // coefficients are given, so beams without aberrations are unaffected.
    if(zernikeCoeffs.empty())
        return 0.0;
    if(beamRadius <= 0.0)
    {
        std::cerr << "Zernike coefficients require a positive 'beamradius'. Exiting." << std::endl;
        exit(1);
    }

    double zernikePhase = 0.0;
    for(const auto& [key, amplitude] : zernikeCoeffs)
    {
        // convert Noll index to the (n,m) pair the polynomials are written in
        std::array<int, 2> zernikeIndices = nollToZernike(key);
        int n = zernikeIndices[0];
        int m = zernikeIndices[1];

        double rho = sqrt(pos[0]*pos[0] + pos[1]*pos[1]) / beamRadius;
        double theta = (rho == 0.0) ? 0.0 : atan2(pos[1], pos[0]);

        zernikePhase += -1.0 * amplitude * 2.0 * pi * Zmn(m, n, rho, theta);
    }
    return zernikePhase;
}

double AISLaserBeam::GetPhi(const doubleThreeVector& pos)
{
    // check the beam type
    if(beamType == "flat_square")
    {
        return phi0;
    }
    else if(beamType == "interpolated")
    {
        // Sampled beam from aisoptics. The exporter writes the slowly varying
        // envelope phase and deliberately excludes the +-kz carrier, which is
        // exactly the convention GetPhi is expected to return.
        doubleThreeVector shiftedPos = ApplyTipTilt(pos);
        double spatiallyVaryingPhase = trilinearInterpolation(
            shiftedPos[0], shiftedPos[1], shiftedPos[2],
            interpolationGrids->x, interpolationGrids->y, interpolationGrids->z,
            interpolationGrids->phase,
            interpolationGrids->x.size(), interpolationGrids->y.size(), interpolationGrids->z.size());

        return phi0 + spatiallyVaryingPhase;
    }
    else if(beamType == "gaussian")
    {
        // define the phase shift for the reflection symmetry
        double reflectionShift;
        if(k[2] > 0)
        {
            reflectionShift = pi;
        }
        else
        {
            reflectionShift = 0.0;
        }

        doubleThreeVector shiftedPos = {pos[0],pos[1],pos[2]};
        if(k[2] < 0)
        {
            shiftedPos[2] = focalLength - pos[2]; // propagating downward
        }
        else if(k[2] > 0)
        {
            shiftedPos[2] = focalLength + pos[2]; // propagating upward
        }

        return phi0 + reflectionShift + gaussianWavefront(shiftedPos, k[2], w0)
                    + GetZernikePhase(pos) + ModePhaseContribution(pos);
    }
    else if(beamType == "confocal")
    {
        // Confocal concave mirror: incoming and outgoing foci coincide at z = focalLength.
        // Both directions share the same local coordinate from the focus, so the intensity
        // profiles are identical. The reflected (+z) beam has the negated spatial phase;
        // no extra pi shift.
        doubleThreeVector shiftedPos = {pos[0], pos[1], focalLength - pos[2]};
        double wavefront = gaussianWavefront(shiftedPos, k[2], w0);
        if(k[2] > 0)
            return phi0 - wavefront + GetZernikePhase(pos) + ModePhaseContribution(pos);
        else
            return phi0 + wavefront + GetZernikePhase(pos) + ModePhaseContribution(pos);
    }
    return phi0; // unreachable: SetBeamType rejects anything else
}

double AISLaserBeam::GetPhiWrapper(double xi, void* params)
{
    getPhiWrapperParams* p = (getPhiWrapperParams*) params;
    
    doubleThreeVector pos = p->pos;
    int index = p->index;
    std::shared_ptr<AISLaserBeam> beam = p->beam;
    // only modify the selected index
    pos[index] = xi;
    return beam->GetPhi(pos);
}

doubleThreeVector AISLaserBeam::GetDelPhi(const doubleThreeVector& pos)
{
    // Local (transverse + longitudinal) correction to the plane-wave
    // wavevector from the beam's curved wavefront (Gouy phase + radius of
    // curvature), i.e. grad(GetPhi(pos)) with the k*z plane-wave term
    // removed (that part is already handled separately via GetK). Mirrors
    // GetPhi's own shiftedPos construction exactly so the two stay consistent.
    if(beamType == "gaussian")
    {
        doubleThreeVector shiftedPos = {pos[0], pos[1], pos[2]};
        double zSign;
        if(k[2] < 0)
        {
            shiftedPos[2] = focalLength - pos[2];
            zSign = -1.0;
        }
        else
        {
            shiftedPos[2] = focalLength + pos[2];
            zSign = 1.0;
        }
        std::array<double,3> grad = gaussianGradientWavefront(shiftedPos, k[2], w0);
        doubleThreeVector modeGrad = ModeGradContribution(pos);
        return {grad[0] + modeGrad[0], grad[1] + modeGrad[1], zSign*grad[2] + modeGrad[2]};
    }
    else if(beamType == "confocal")
    {
        doubleThreeVector shiftedPos = {pos[0], pos[1], focalLength - pos[2]};
        std::array<double,3> grad = gaussianGradientWavefront(shiftedPos, k[2], w0);
        double sign = (k[2] > 0) ? -1.0 : 1.0; // matches GetPhi's +-wavefront choice
        doubleThreeVector modeGrad = ModeGradContribution(pos);
        return {sign*grad[0] + modeGrad[0], sign*grad[1] + modeGrad[1],
                -sign*grad[2] + modeGrad[2]};
    }
    // KNOWN ISSUE (v0.0.2): beamType == "interpolated" returns a zero gradient,
    // so sampled beams carry no wavefront-gradient recoil. GetPhiWrapper and
    // gsl_deriv.h are retained for the numerical implementation; note that
    // trilinear interpolation is only C0, so the differentiated gradient is
    // discontinuous across cell boundaries and a smoother interpolant may be
    // needed before this is trustworthy. Zernike aberrations have the same gap:
    // GetPhi applies them, GetDelPhi does not. See KNOWN_ISSUES.md.
    return {0., 0., 0.};
}

double AISLaserBeam::GetRabiFreq(const doubleThreeVector& pos, const __float128& t0, const __float128& t)
{
    // check the beam type
    double effectiveRabiFreq = 0.0;
    if(beamType == "flat_square")
    {
        effectiveRabiFreq = rabiFreq; // constant rabi freq
    }
    else if(beamType == "interpolated")
    {
        // aisoptics writes field amplitude relative to the configured Rabi
        // frequency, so the sampled value multiplies rabiFreq directly.
        doubleThreeVector shiftedPos = ApplyTipTilt(pos);
        double amplitude = trilinearInterpolation(
            shiftedPos[0], shiftedPos[1], shiftedPos[2],
            interpolationGrids->x, interpolationGrids->y, interpolationGrids->z,
            interpolationGrids->amplitude,
            interpolationGrids->x.size(), interpolationGrids->y.size(), interpolationGrids->z.size());

        effectiveRabiFreq = amplitude * rabiFreq;
    }
    else if(beamType == "gaussian")
    {
        double zPos;
        if (k[2] > 0)
        {
            zPos = focalLength + pos[2];
        }
        else if (k[2] < 0)
        {
            zPos = focalLength - pos[2];
        }
        doubleThreeVector shiftedPos = {pos[0], pos[1], zPos};
        effectiveRabiFreq = rabiFreq * gaussianEnvelope(shiftedPos, w0, t0, t)
                                     * ModeAmplitudeContribution(pos);
    }
    else if(beamType == "confocal")
    {
        // Same intensity profile for both directions: both measure distance from the shared focus.
        doubleThreeVector shiftedPos = {pos[0], pos[1], focalLength - pos[2]};
        effectiveRabiFreq = rabiFreq * gaussianEnvelope(shiftedPos, w0, t0, t)
                                     * ModeAmplitudeContribution(pos);
    }
    return effectiveRabiFreq;
}

double AISLaserBeam::GetW0()
{
    return w0;
}
void AISLaserBeam::SetW0(double w0Value)
{
    w0 = w0Value;
}

double AISLaserBeam::GetZLaser()
{
    return zLaser;
}
void AISLaserBeam::SetZLaser(double zLaserValue)
{
    zLaser = zLaserValue;
}

double AISLaserBeam::GetFocalLength()
{
    return focalLength;
}
void AISLaserBeam::SetFocalLength(double focalLengthValue)
{
    focalLength = focalLengthValue;
}

void AISLaserBeam::SetZernikeCoeffs(std::map<int, double> zernikeCoeffsValue)
{
    zernikeCoeffs = zernikeCoeffsValue;
}
std::map<int, double> AISLaserBeam::GetZernikeCoeffs()
{
    return zernikeCoeffs;
}

void AISLaserBeam::SetBeamRadius(double beamRadiusValue)
{
    beamRadius = beamRadiusValue;
}
double AISLaserBeam::GetBeamRadius()
{
    return beamRadius;
}

void AISLaserBeam::SetBeamType(std::string beamTypeValue)
{
    // check if 'flat_square', 'gaussian', 'confocal' or 'interpolated'
    if(beamTypeValue == "flat_square" || beamTypeValue == "gaussian" ||
       beamTypeValue == "confocal"    || beamTypeValue == "interpolated")
    {
        beamType = beamTypeValue;
    }
    else
    {
        std::cerr << "Beam type " << beamTypeValue << " not recognized. Exiting." << std::endl;
        exit(1);
    }
}
std::string AISLaserBeam::GetBeamType()
{
    return beamType;
}

void AISLaserBeam::SetBaselineLength(double baselineLengthValue)
{
    baselineLength = baselineLengthValue;
}
double AISLaserBeam::GetBaselineLength()
{
    return baselineLength;
}

void AISLaserBeam::SetTipTilt(double tiptiltXValue, double tiptiltYValue)
{
    tiptiltX = tiptiltXValue;
    tiptiltY = tiptiltYValue;
}

std::pair<double, double> AISLaserBeam::GetTipTilt()
{
    return {tiptiltX, tiptiltY};
}

bool AISLaserBeam::HasInterpolationGrids() const
{
    return interpolationGrids != nullptr;
}

std::map<std::string, std::shared_ptr<const AISLaserBeam::InterpolationGrids>>&
AISLaserBeam::InterpolationCache()
{
    static std::map<std::string, std::shared_ptr<const InterpolationGrids>> cache;
    return cache;
}

size_t AISLaserBeam::InterpolationCacheSize()
{
    return InterpolationCache().size();
}

void AISLaserBeam::SetInterpolationGrids(std::string filename)
{
    if(filename.empty())
    {
        std::cerr << "Beam type 'interpolated' requires a file in "
                     "'beaminterpolationparamsfilenames'. Exiting." << std::endl;
        exit(1);
    }

    // A real pulse sequence typically has hundreds to thousands of pulses
    // sharing one or two beam files; load each distinct file once and share
    // it by shared_ptr across every beam that references it; see the
    // InterpolationGrids/InterpolationCache comment in the header for why
    // a naive per-beam copy is not just wasteful but can exhaust memory.
    auto& cache = InterpolationCache();
    auto cached = cache.find(filename);
    if(cached != cache.end())
    {
        interpolationGrids = cached->second;
        return;
    }

    auto grids = std::make_shared<InterpolationGrids>();

    try
    {
        H5::H5File file(filename, H5F_ACC_RDONLY);

        // Refuse a layout this reader was not written against. aisoptics stamps
        // export_format on every file it writes; a mismatch here is what stops a
        // future layout change from silently producing wrong phases.
        if(file.attrExists("export_format"))
        {
            H5::Attribute attr = file.openAttribute("export_format");
            std::string exportFormat;
            attr.read(attr.getStrType(), exportFormat);
            if(exportFormat != "aispp_current_interpolation_hdf5")
            {
                std::cerr << "Beam file " << filename << " has export_format '" << exportFormat
                          << "', expected 'aispp_current_interpolation_hdf5'. Exiting." << std::endl;
                exit(1);
            }
        }
        else
        {
            std::cerr << "Warning: beam file " << filename << " has no 'export_format' attribute; "
                      << "assuming the aisoptics total-field layout." << std::endl;
        }

        auto readVector = [&file](const std::string& name, std::vector<double>& out)
        {
            H5::DataSet dataset = file.openDataSet(name);
            H5::DataSpace dataspace = dataset.getSpace();
            hsize_t dims[1];
            dataspace.getSimpleExtentDims(dims);
            out.resize(dims[0]);
            dataset.read(out.data(), H5::PredType::NATIVE_DOUBLE);
        };

        readVector("x", grids->x);
        readVector("y", grids->y);
        readVector("z", grids->z);
        readVector("phase", grids->phase);
        readVector("amplitude", grids->amplitude);

        file.close();
    }
    catch(const H5::Exception& e)
    {
        std::cerr << "Failed to read beam file " << filename << ": "
                  << e.getDetailMsg() << ". Exiting." << std::endl;
        exit(1);
    }

    // trilinearInterpolation needs at least one cell along each axis, and the flat
    // value arrays must match the product of the axis lengths under the C-order
    // index i*Ny*Nz + j*Nz + k that aisoptics documents.
    if(grids->x.size() < 2 || grids->y.size() < 2 || grids->z.size() < 2)
    {
        std::cerr << "Beam file " << filename << " needs at least 2 points on each axis. Exiting." << std::endl;
        exit(1);
    }

    const size_t expected = grids->x.size() * grids->y.size() * grids->z.size();
    if(grids->phase.size() != expected || grids->amplitude.size() != expected)
    {
        std::cerr << "Beam file " << filename << ": phase/amplitude length ("
                  << grids->phase.size() << "/" << grids->amplitude.size()
                  << ") does not match Nx*Ny*Nz = " << expected << ". Exiting." << std::endl;
        exit(1);
    }

    cache[filename] = grids;
    interpolationGrids = grids;
}



void AISLaserBeam::SetFourierMode(const FourierMode& mode)
{
    fourierMode = mode;
}

FourierMode AISLaserBeam::GetFourierMode() const
{
    return fourierMode;
}

void AISLaserBeam::SetWavefrontGradientEnabled(bool enabled)
{
    wavefrontGradientEnabled = enabled;
}

bool AISLaserBeam::GetWavefrontGradientEnabled() const
{
    return wavefrontGradientEnabled;
}

void AISLaserBeam::CarrierCoordinates(const doubleThreeVector& pos,
                                      doubleThreeVector& shiftedPos,
                                      double& dzFromMirror,
                                      double& zSign) const
{
    // Reproduces GetPhi's own shiftedPos construction for each beam type, so
    // the mode is evaluated in exactly the coordinate the carrier uses.
    if(beamType == "confocal")
    {
        shiftedPos = {pos[0], pos[1], focalLength - pos[2]};
        zSign = -1.0;
    }
    else if(k[2] < 0)
    {
        shiftedPos = {pos[0], pos[1], focalLength - pos[2]};
        zSign = -1.0;
    }
    else
    {
        shiftedPos = {pos[0], pos[1], focalLength + pos[2]};
        zSign = 1.0;
    }
    dzFromMirror = pos[2] - fourierMode.mirrorZ;
}

void AISLaserBeam::FourierModeContributions(const doubleThreeVector& pos,
                                            double& phase, double& amplitudeFactor) const
{
    phase = 0.0;
    amplitudeFactor = 1.0;
    if(!fourierMode.active())
    {
        return;
    }
    doubleThreeVector shiftedPos;
    double dzFromMirror, zSign;
    CarrierCoordinates(pos, shiftedPos, dzFromMirror, zSign);
    // confocal with kz > 0 represents the outgoing beam with a negated spatial
    // phase, i.e. the conjugate carrier; building sidebands from the
    // unconjugated one flips the sign of the curvature phase difference across
    // the walk-off and costs several percent for a tightly focused beam.
    const bool conjugate = (beamType == "confocal" && k[2] > 0);
    fourierModeFirstOrder(shiftedPos, pos[0], pos[1], dzFromMirror, k[2], w0,
                          fourierMode, conjugate, phase, amplitudeFactor);
}

double AISLaserBeam::GetFourierModePhase(const doubleThreeVector& pos) const
{
    double phase, amplitudeFactor;
    FourierModeContributions(pos, phase, amplitudeFactor);
    return phase;
}

double AISLaserBeam::GetFourierModeAmplitudeFactor(const doubleThreeVector& pos) const
{
    double phase, amplitudeFactor;
    FourierModeContributions(pos, phase, amplitudeFactor);
    return amplitudeFactor;
}

doubleThreeVector AISLaserBeam::GetFourierModeGradPhase(const doubleThreeVector& pos) const
{
    // Gradient of the phase the mode adds, by central differences on an
    // analytic, infinitely differentiable function. This is the crucial
    // difference from the sampled-beam case: trilinear interpolation is only
    // C0, so differentiating it gives a gradient that jumps across cell
    // boundaries (KNOWN_ISSUES.md). Here the underlying function is smooth, so
    // a centred difference is well conditioned.
    if(!fourierMode.active() || !wavefrontGradientEnabled)
    {
        return {0.0, 0.0, 0.0};
    }
    const double q = fourierMode.q();
    const double h = (q > 0.0) ? std::min(1.0e-6, 0.001 * 2.0 * pi / q) : 1.0e-6;

    doubleThreeVector grad{0.0, 0.0, 0.0};
    for(int axis = 0; axis < 3; ++axis)
    {
        doubleThreeVector plus = pos, minus = pos;
        plus[axis]  += h;
        minus[axis] -= h;
        grad[axis] = (GetFourierModePhase(plus) - GetFourierModePhase(minus)) / (2.0 * h);
    }
    return grad;
}


double AISLaserBeam::ModePhaseContribution(const doubleThreeVector& pos) const
{
    // The ripple is a property of the MIRROR, imprinted on reflection, so it
    // rides only on the outgoing (+z) beam. The downgoing beam is unaberrated,
    // which is the same split the (up_file, down_file) pair expresses for
    // sampled beams.
    if(!fourierMode.active() || k[2] <= 0)
    {
        return 0.0;
    }
    // The ripple's sign is a property of the mirror surface, not of how the
    // carrier's Gouy/curvature phase is bookkept, so it is NOT negated for
    // confocal. Validated against the sampled-beam route, which reproduces the
    // reference project's PSMAPs: negating it here flipped the sign of the
    // wavefront phase relative to that reference.
    return GetFourierModePhase(pos);
}

double AISLaserBeam::ModeAmplitudeContribution(const doubleThreeVector& pos) const
{
    if(!fourierMode.active() || k[2] <= 0)
    {
        return 1.0;
    }
    return GetFourierModeAmplitudeFactor(pos);
}


doubleThreeVector AISLaserBeam::ModeGradContribution(const doubleThreeVector& pos) const
{
    // Gradient of exactly what ModePhaseContribution adds, so GetDelPhi stays
    // the gradient of GetPhi. This is the term GetDelPhi has always returned as
    // zero for aberrated beams; with an analytic mode it is available, and
    // SetWavefrontGradientEnabled(false) restores the old behaviour so runs
    // with and without the gradient recoil can be compared.
    if(!fourierMode.active() || k[2] <= 0 || !wavefrontGradientEnabled)
    {
        return {0.0, 0.0, 0.0};
    }
    return GetFourierModeGradPhase(pos);
}
