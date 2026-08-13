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
            interpolationGridX, interpolationGridY, interpolationGridZ,
            phaseInterpolationGridValues,
            interpolationGridX.size(), interpolationGridY.size(), interpolationGridZ.size());

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
                    + GetZernikePhase(pos);
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
            return phi0 - wavefront + GetZernikePhase(pos);
        else
            return phi0 + wavefront + GetZernikePhase(pos);
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
        return {grad[0], grad[1], zSign * grad[2]};
    }
    else if(beamType == "confocal")
    {
        doubleThreeVector shiftedPos = {pos[0], pos[1], focalLength - pos[2]};
        std::array<double,3> grad = gaussianGradientWavefront(shiftedPos, k[2], w0);
        double sign = (k[2] > 0) ? -1.0 : 1.0; // matches GetPhi's +-wavefront choice
        return {sign*grad[0], sign*grad[1], -sign*grad[2]};
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
            interpolationGridX, interpolationGridY, interpolationGridZ,
            amplitudeInterpolationGridValues,
            interpolationGridX.size(), interpolationGridY.size(), interpolationGridZ.size());

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
        effectiveRabiFreq = rabiFreq * gaussianEnvelope(shiftedPos, w0, t0, t);
    }
    else if(beamType == "confocal")
    {
        // Same intensity profile for both directions: both measure distance from the shared focus.
        doubleThreeVector shiftedPos = {pos[0], pos[1], focalLength - pos[2]};
        effectiveRabiFreq = rabiFreq * gaussianEnvelope(shiftedPos, w0, t0, t);
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
    return !interpolationGridX.empty();
}

void AISLaserBeam::SetInterpolationGrids(std::string filename)
{
    if(filename.empty())
    {
        std::cerr << "Beam type 'interpolated' requires a file in "
                     "'beaminterpolationparamsfilenames'. Exiting." << std::endl;
        exit(1);
    }

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

        readVector("x", interpolationGridX);
        readVector("y", interpolationGridY);
        readVector("z", interpolationGridZ);
        readVector("phase", phaseInterpolationGridValues);
        readVector("amplitude", amplitudeInterpolationGridValues);

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
    if(interpolationGridX.size() < 2 || interpolationGridY.size() < 2 || interpolationGridZ.size() < 2)
    {
        std::cerr << "Beam file " << filename << " needs at least 2 points on each axis. Exiting." << std::endl;
        exit(1);
    }

    const size_t expected = interpolationGridX.size() * interpolationGridY.size() * interpolationGridZ.size();
    if(phaseInterpolationGridValues.size() != expected || amplitudeInterpolationGridValues.size() != expected)
    {
        std::cerr << "Beam file " << filename << ": phase/amplitude length ("
                  << phaseInterpolationGridValues.size() << "/" << amplitudeInterpolationGridValues.size()
                  << ") does not match Nx*Ny*Nz = " << expected << ". Exiting." << std::endl;
        exit(1);
    }
}

