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

doubleThreeVector AISLaserBeam::GetK()
{
    return this->k;
}

void AISLaserBeam::SetK(doubleThreeVector k)
{
    this->k = k;
}

__float128 AISLaserBeam::GetOmega()
{
    return this->omega;
}

void AISLaserBeam::SetOmega(__float128 omega)
{
    this->omega=omega;
}

double AISLaserBeam::GetPhi(const doubleThreeVector& pos)
{
    return phi0 + (*wavefrontFunction)(pos);
}

doubleThreeVector AISLaserBeam::GetDelPhi(const doubleThreeVector& pos)
{
    // Check if delWavefrontFunction is null
    if (delWavefrontFunction == nullptr) {
        std::cerr << "Error: delWavefrontFunction is not initialized!" << std::endl;
        throw std::runtime_error("delWavefrontFunction is not initialized");
    }
    return (*delWavefrontFunction)(pos);
}

double AISLaserBeam::GetRabiFreq(const doubleThreeVector& pos, const __float128& t0, const __float128& t)
{
    return this->rabiFreq * (*rabifreqFunction)(pos, t0, t); // central rabi freq times envelope
}