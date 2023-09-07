#ifndef AISLINEARGRAVITYPROPAGATOR_HH
#define AISLINEARGRAVITYPROPAGATOR_HH

#include "AISAbcdPropagator.hh"
#include "AISConstants.hh"

class AISLinearGravityPropagator : public AISAbcdPropagator 
{
public :
    AISLinearGravityPropagator();
    ~AISLinearGravityPropagator();

protected:
    quad3x3Matrix alpha = {};
    quad3x3Matrix gamma = {{ {{-g/radiusEarth,      0.0q,            0.0q     }},
                             {{     0.0q     , -g/radiusEarth,       0.0q     }},
                             {{     0.0q     ,      0.0q,      2*g/radiusEarth}} }};

    quadThreeVector gVector = {0.0q, 0.0q, -g};
}

AISLinearGravityPropagator::AISLinearGravityPropagator(){}
AISLinearGravityPropagator::~AISLinearGravityPropagator(){}

#endif