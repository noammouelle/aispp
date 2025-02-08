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
        // compute the gaussian phase
        double zShift;
        if (k[2] > 0)
        {
            zShift = focalLength + zLaser;
        }
        else
        {
            zShift = zLaser - focalLength + baselineLength;
        }
        doubleThreeVector shiftedPos = {pos[0], pos[1], pos[2] - zShift}; // shift the position to account for position of the lens
                                                                          // and the laser focal length

        double gaussianPhase = gaussianWavefront(shiftedPos, k[2], w0);
        
        // compute the zernike polynomial phase
        double zernikePhase = 0.0;
        for(const auto& [key, amplitude] : zernikeCoeffs)
        {
            // convert OSA/ANSI Zernike polynomial index to n,m
            std::array<int, 2> zernikeIndices = nollToZernike(key);
            int n = zernikeIndices[0];
            int m = zernikeIndices[1];
            // compute rho, theta
            double rho = sqrt(pos[0]*pos[0] + pos[1]*pos[1]) / beamRadius;
            double theta = atan2(pos[1],pos[0]) + pi;//atan2(pos[1], pos[0]);
            // compute the zernike polynomial phase (same convention as in https://opticspy.github.io/lightpipes/command-reference.html#LightPipes.Zernike)
            double prefactor = sqrt((2.0 * n + 2.0) / (1.0 + (m == 0)));

            // define the sign (if propagating in the negative z direction, the sign is inverted)
            double sign;
            if(k[2] < 0)
            {
                sign = -1.0;
            }
            else
            {
                sign = 1.0;
            }

            zernikePhase +=  -1.0 * sign * prefactor * amplitude * 2.0 * pi * Zmn(n, m, rho, theta);
        }

        return phi0 + zernikePhase;
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
    // check the beam type
    if(beamType == "flat_square")
    {
        return {0.0, 0.0, 0.0};
    }
    else if(beamType == "gaussian")
    {
        doubleThreeVector delPhi = {0.0, 0.0, 0.0};
        double result, abserr;

        for (int i = 0; i < 2; ++i)
        {
            getPhiWrapperParams params;
            params.pos = pos;
            params.index = i;
            params.beam = std::make_shared<AISLaserBeam>(*this);

            gsl_function F;
            F.params = &params;

            F.function = GetPhiWrapper;
            gsl_deriv_central(&F, pos[i], 1e-8, &result, &abserr);
            delPhi[i] = result;
        }

        return delPhi;
    }
}

double AISLaserBeam::GetRabiFreq(const doubleThreeVector& pos, const __float128& t0, const __float128& t)
{
    // check the beam type
    if(beamType == "flat_square")
    {
        return rabiFreq;
    }
    else if(beamType == "gaussian")
    {
        double zShift;
        if (k[2] > 0)
        {
            zShift = focalLength + zLaser;
        }
        else
        {
            zShift = zLaser - focalLength + baselineLength;
        }
        doubleThreeVector shiftedPos = {pos[0], pos[1], pos[2] - zShift}; // shift the position to account for position of the lens
                                                                          // and the laser focal length
        double effectiveRabiFreq = rabiFreq * gaussianEnvelope(shiftedPos, w0, t0, t); // central rabi freq times envelope

        return effectiveRabiFreq; // central rabi freq times envelope
    }
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
    if(beamTypeValue == "flat_square" || beamTypeValue == "gaussian")
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

