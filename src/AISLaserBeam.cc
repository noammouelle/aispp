#include "AISLaserBeam.cc"

AISLaserBeam::AISLaserBeam(doubleThreeVector k, __float128 omega)
{
    this->k = k
    this->omega = omega
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