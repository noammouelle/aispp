#ifndef AISLINEARGRAVITYPROPAGATOR_HH
#define AISLINEARGRAVITYPROPAGATOR_HH

#include "AISAbcdPropagator.hh"
#include "AISConstants.hh"

class AISLinearGravityPropagator : public AISAbcdPropagator 
{
public :
    AISLinearGravityPropagator(__float128 dt);
    ~AISLinearGravityPropagator();
};

#endif