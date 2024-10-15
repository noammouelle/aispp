#include "AISLinearGravityPropagator.hh"

AISLinearGravityPropagator::AISLinearGravityPropagator()
{
    // define here the gravity gradient tensor
    GGtensor = {{{g / radiusEarth, 0.0, 0.0},
                                {0.0, g / radiusEarth, 0.0},
                                {0.0, 0.0, - 2.0 * g / radiusEarth}}};
}

AISLinearGravityPropagator::~AISLinearGravityPropagator()
{}

double AISLinearGravityPropagator::get_U(const doubleThreeVector& pos, const doubleThreeVector& vel)
{
    //return massSr87 * g * pos[2] + 0.5 * massSr87 * dotProduct(dotProduct(GGtensor, pos), pos);
    return g * pos[2] + 0.5 * dotProduct(dotProduct(GGtensor, pos), pos);
}

doubleThreeVector AISLinearGravityPropagator::get_dUdx(const doubleThreeVector& pos, const doubleThreeVector& vel)
{
    doubleThreeVector term1 = {0.0, 0.0, massSr87 * g};
    doubleThreeVector term2 = scalarMultiply(dotProduct(GGtensor, pos), massSr87);
    return matrixAdd(term1, term2);
}

doubleThreeVector AISLinearGravityPropagator::get_dUdp(const doubleThreeVector& pos, const doubleThreeVector& vel)
{
    return {0.0, 0.0, 0.0};
}