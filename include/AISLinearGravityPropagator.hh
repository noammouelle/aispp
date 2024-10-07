#ifndef AISLINEARGRAVITYPROPAGATOR_HH
#define AISLINEARGRAVITYPROPAGATOR_HH

#include "AISKinematicPropagator.hh"
#include "AISConstants.hh"
#include "AISUtilities.hh"

class AISLinearGravityPropagator : public AISKinematicPropagator
{
public :
    AISLinearGravityPropagator();
    ~AISLinearGravityPropagator();
};

#endif