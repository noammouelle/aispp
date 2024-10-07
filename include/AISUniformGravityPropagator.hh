#ifndef AISUNIFORMGRAVITYPROPAGATOR_HH
#define AISUNIFORMGRAVITYPROPAGATOR_HH

#include "AISKinematicPropagator.hh"
#include "AISConstants.hh"

class AISUniformGravityPropagator : public AISKinematicPropagator
{
public :
    AISUniformGravityPropagator();
    ~AISUniformGravityPropagator();
};

#endif