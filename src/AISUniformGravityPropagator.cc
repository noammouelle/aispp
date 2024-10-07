#include "AISUniformGravityPropagator.hh"

AISUniformGravityPropagator::AISUniformGravityPropagator()
{
}

AISUniformGravityPropagator::~AISUniformGravityPropagator()
{}

double AISUniformGravityPropagator::get_U(const doubleThreeVector& pos, const doubleThreeVector& vel)
{
    return massSr87 * g * pos[2];
}

doubleThreeVector AISUniformGravityPropagator::get_dUdx(const doubleThreeVector& pos, const doubleThreeVector& vel)
{
    return {0.0, 0.0, massSr87 * g};
}

doubleThreeVector AISUniformGravityPropagator::get_dUdp(const doubleThreeVector& pos, const doubleThreeVector& vel)
{
    return {0.0, 0.0, 0.0};
}