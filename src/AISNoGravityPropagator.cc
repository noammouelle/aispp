#include "AISNoGravityPropagator.hh"

AISNoGravityPropagator::AISNoGravityPropagator()
{
}

AISNoGravityPropagator::~AISNoGravityPropagator()
{}

double AISNoGravityPropagator::get_U(const doubleThreeVector& pos, const doubleThreeVector& vel)
{
    return 0;
}

doubleThreeVector AISNoGravityPropagator::get_dUdx(const doubleThreeVector& pos, const doubleThreeVector& vel)
{
    return {0.0, 0.0, 0.0};
}

doubleThreeVector AISNoGravityPropagator::get_dUdp(const doubleThreeVector& pos, const doubleThreeVector& vel)
{
    return {0.0, 0.0, 0.0};
}