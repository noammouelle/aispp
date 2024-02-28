#ifndef AISATOMINTERFEROMETERPARAMS_HH
#define AISATOMINTERFEROMETERPARAMS_HH

#include "AISUtilities.hh"

class AISAtomInterferometerParams
{
public:
    AISAtomInterferometerParams();
    ~AISAtomInterferometerParams();

    /* CLOUD PARAMETERS */
    doubleThreeVector initialPosition;
    doubleThreeVector initialVelocity;
    double cloudTemperature;
    double cloudWidth;
    int nAtoms;

    /* BEAM PARAMETERS */
    double beamRadius;
    double laserPhase;
    double rabiFrequency;
    doubleThreeVector psrGradient;
    doubleThreeVector aberrationK;
    double aberrationAmplitude;

    /* SEQUENCE PARAMETERS */
    __float128 initialPropagationTime;
    __float128 finalPropagationTime;
    __float128 interogationTime;
    int nSteps;

    /* LMT PARAMS */
    int lmtOrder;
    double lmtDelayTime;

    /* DETECTOR PARAMETERS */
    double coherenceLength;
};

#endif