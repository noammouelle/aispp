#ifndef AISLINEARGRAVITYPROPAGATOR_HH
#define AISLINEARGRAVITYPROPAGATOR_HH

#include "AISFreePropagator.hh"
#include "AISConstants.hh"

class AISUniformGravityPropagator : public AISKinematicPropagator
{
public :
    AISUniformGravityPropagator();
    ~AISUniformGravityPropagator();
};

#endif