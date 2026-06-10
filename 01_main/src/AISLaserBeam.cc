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

double AISLaserBeam::GetPhi(const doubleThreeVector& pos)
{
    // check the beam type
    if(beamType == "flat_square")
    {
        return phi0;
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

        return phi0 + reflectionShift + gaussianWavefront(shiftedPos, k[2], w0);
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
            return phi0 - wavefront;
        else
            return phi0 + wavefront;
    }
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
    return {0., 0., 0.}; // TODO: implement 
}

double AISLaserBeam::GetRabiFreq(const doubleThreeVector& pos, const __float128& t0, const __float128& t)
{
    // check the beam type
    double effectiveRabiFreq = 0.0;
    if(beamType == "flat_square")
    {
        effectiveRabiFreq = rabiFreq; // constant rabi freq
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
    // check if 'flat_square' or 'gaussian'
    if(beamTypeValue == "flat_square" || beamTypeValue == "gaussian" || beamTypeValue == "confocal")
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

