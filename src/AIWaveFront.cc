#include "AISWaveFront.hh"
#include "AISUtilities.hh"
#include "AISConstants.hh"

AISWaveFront::AISWaveFront(double (*aberrationFunction)(const doubleThreeVector&), 
                           double (*aberrationFunctionDerivative)(const doubleThreeVector&),
                           doubleThreeVector psrGradient, double laserPhase)
{
    fAberrationFunction = aberrationFunction;
    fAberrationFunctionDerivative = aberrationFunctionDerivative;
    fPsrGradient = psrGradient;
    fLaserPhase = laserPhase;
}

AISWaveFront::~AISWaveFront(){}

double AISWaveFront::GetValue(const  doubleThreeVector& pos)
{
    // compute the linear component
    double linearComponent = dotProduct(fPsrGradient, pos) + fLaserPhase;
    // compute the aberration
    double aberration = fAberrationFunction(pos);
    // return the wavefront value
    return linearComponent + aberration;
}

// define a few aberrationFunctions
double zeroAberrationFunction(const doubleThreeVector& pos)
{
    return 0.;
}
double sineAberrationFunction(const doubleThreeVector& pos, double A, doubleThreeVector k)
{
    return A * sin(dotProduct(k, pos));
}