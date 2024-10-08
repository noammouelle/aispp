#ifndef AISLINEARGRAVITYPROPAGATOR_HH
#define AISLINEARGRAVITYPROPAGATOR_HH

#include "AISKinematicPropagator.hh"
#include "AISConstants.hh"
#include "AISUtilities.hh"

class AISLinearGravityPropagator : public AISKinematicPropagator
{
protected:
    double get_U(const doubleThreeVector& pos, const doubleThreeVector& vel) override;
    doubleThreeVector get_dUdx(const doubleThreeVector& pos, const doubleThreeVector& vel) override;
    doubleThreeVector get_dUdp(const doubleThreeVector& pos, const doubleThreeVector& vel) override;

    double3x3Matrix GGtensor;
    
public :
    AISLinearGravityPropagator();
    ~AISLinearGravityPropagator();
};

#endif