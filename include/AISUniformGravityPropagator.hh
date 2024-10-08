#ifndef AISUNIFORMGRAVITYPROPAGATOR_HH
#define AISUNIFORMGRAVITYPROPAGATOR_HH

#include "AISKinematicPropagator.hh"
#include "AISConstants.hh"

class AISUniformGravityPropagator : public AISKinematicPropagator
{
protected:
    double get_U(const doubleThreeVector& pos, const doubleThreeVector& vel) override;
    doubleThreeVector get_dUdx(const doubleThreeVector& pos, const doubleThreeVector& vel) override;
    doubleThreeVector get_dUdp(const doubleThreeVector& pos, const doubleThreeVector& vel) override;
public :
    AISUniformGravityPropagator();
    ~AISUniformGravityPropagator();
};

#endif