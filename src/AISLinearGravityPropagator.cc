#include "AISLinearGravityPropagator.hh"
#include "AISUtilities.hh"

#include <iostream>

AISLinearGravityPropagator::AISLinearGravityPropagator(__float128 dt) : AISAbcdPropagator(dt)
{
    deltaTime   = dt;
    deltaTime64 = convertScalarToDouble(dt);
    //deltaTime64 = dt;

    alpha = {{ {{0.0, 0.0, 0.0}},
               {{0.0, 0.0, 0.0}},
               {{0.0, 0.0, 0.0}} }};

    gamma = {{ {{-g/radiusEarth,      0.0,            0.0     }},
               {{     0.0     , -g/radiusEarth,       0.0     }},
               {{     0.0     ,      0.0,      2*g/radiusEarth}} }};

    gVector = {0.0, 0.0, -g};

    setABCDXiPhi();
}

AISLinearGravityPropagator::~AISLinearGravityPropagator(){}