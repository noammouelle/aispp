#include "AISLaserBeam.hh"

AISLaserBeam::AISLaserBeam(doubleThreeVector k, __float128 omega,
                           wavefrontFunctionType wavefrontFunction,
                           delWavefrontFunctionType delWavefrontFunction,
                           rabifreqFunctionType rabiFreqFunction)
{
    this->k = k;
    this->omega = omega;
    this->wavefrontFunction = wavefrontFunction;
    this->delWavefrontFunction = delWavefrontFunction;
    this->rabifreqFunction = rabiFreqFunction;

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

double AISLaserBeam::GetPhi(doubleThreeVector pos)
{
    return wavefrontFunction(pos);
}

doubleThreeVector AISLaserBeam::GetDelPhi(doubleThreeVector pos)
{
    return delWavefrontFunction(pos);
}

double AISLaserBeam::GetRabiFreq(doubleThreeVector pos, __float128 t)
{
    return rabifreqFunction(pos, t);
}