#ifndef AISNOGRAVITYPROPAGATOR_HH
#define AISNOGRAVITYPROPAGATOR_HH

#include "AISKinematicPropagator.hh"
#include "AISConstants.hh"

class AISNoGravityPropagator : public AISKinematicPropagator
{
protected:
    double get_U(const doubleThreeVector& pos, const doubleThreeVector& vel) override;
    doubleThreeVector get_dUdx(const doubleThreeVector& pos, const doubleThreeVector& vel) override;
    doubleThreeVector get_dUdp(const doubleThreeVector& pos, const doubleThreeVector& vel) override;
public :
    AISNoGravityPropagator();
    ~AISNoGravityPropagator();
};

#endif