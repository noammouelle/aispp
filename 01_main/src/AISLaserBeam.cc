#include "AISLaserBeam.hh"

AISLaserBeam::AISLaserBeam(doubleThreeVector aK, __float128 aOmega, double aRabiFreq, double phi0,
                           std::shared_ptr<wavefrontFunctionType> aWavefrontFunction,
                           std::shared_ptr<delWavefrontFunctionType> aDelWavefrontFunction,
                           std::shared_ptr<rabifreqFunctionType> aRabiFreqFunction) : k(aK), omega(aOmega), rabiFreq(aRabiFreq), phi0(phi0),
                            wavefrontFunction(aWavefrontFunction), delWavefrontFunction(aDelWavefrontFunction),
                            rabifreqFunction(aRabiFreqFunction)
{
}

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
    double zShift;
    if (k[2] > 0)
    {
        zShift = focalLength + zLaser;
    }
    else
    {
        zShift = - focalLength + zLaser;
    }
    doubleThreeVector shiftedPos = {pos[0], pos[1], pos[2] - zShift}; // shift the position to account for position of the lens
                                                                      // and the laser focal length
    return phi0 + (*wavefrontFunction)(shiftedPos, k[2], w0);
}

doubleThreeVector AISLaserBeam::GetDelPhi(const doubleThreeVector& pos)
{
    // Check if delWavefrontFunction is null
    if (delWavefrontFunction == nullptr) {
        std::cerr << "Error: delWavefrontFunction is not initialized!" << std::endl;
        throw std::runtime_error("delWavefrontFunction is not initialized");
    }
    double zShift;
    if (k[2] > 0)
    {
        zShift = focalLength + zLaser;
    }
    else
    {
        zShift = - focalLength + zLaser;
    }
    doubleThreeVector shiftedPos = {pos[0], pos[1], pos[2] - zShift}; // shift the position to account for position of the lens
                                                                      // and the laser focal length
    return (*delWavefrontFunction)(shiftedPos, k[2], w0);
}

double AISLaserBeam::GetRabiFreq(const doubleThreeVector& pos, const __float128& t0, const __float128& t)
{
    double zShift;
    if (k[2] > 0)
    {
        zShift = focalLength + zLaser;
    }
    else
    {
        zShift = - focalLength + zLaser;
    }
    doubleThreeVector shiftedPos = {pos[0], pos[1], pos[2] - zShift}; // shift the position to account for position of the lens
                                                                    // and the laser focal length
    double effectiveRabiFreq = rabiFreq * (*rabifreqFunction)(shiftedPos, w0, t0, t); // central rabi freq times envelope

    return effectiveRabiFreq; // central rabi freq times envelope
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