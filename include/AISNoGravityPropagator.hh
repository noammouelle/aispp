#ifndef AISNOGRAVITYPROPAGATOR_HH
#define AISNOGRAVITYPROPAGATOR_HH

#include "AISKinematicPropagator.hh"
#include "AISConstants.hh"

class AISNoGravityPropagator : public AISKinematicPropagator
{
public :
    AISNoGravityPropagator();
    ~AISNoGravityPropagator();
};

#endif