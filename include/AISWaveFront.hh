#ifndef AISWAVEFRONT_HH
#define AISWAVEFRONT_HH

#include <array>
#include "AISUtilities.hh"

class AISWaveFront
{
private:
    double (*fAberrationFunction)(const doubleThreeVector&);
    double (*fAberrationFunctionDerivative)(const doubleThreeVector&);
    doubleThreeVector fPsrGradient;
    double fLaserPhase;
public:
    AISWaveFront(double (*aberrationFunction)(const doubleThreeVector&), 
                 double (*aberrationFunctionDerivative)(const doubleThreeVector&),
                 doubleThreeVector psrGradient, double laserPhase);
    ~AISWaveFront();

    double GetValue(const doubleThreeVector& pos);
};

// define a few aberrationFunctions
double zeroAberrationFunction(const doubleThreeVector& pos);
double sineAberrationFunction(const doubleThreeVector& pos, double A, double k);

#endif