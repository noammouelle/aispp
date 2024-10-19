#ifndef AISLASERBEAM_HH
#define AISLASERBEAM_HH

#include "AISConstants.hh"
#include "AISUtilities.hh"

class AISLaserBeam
{
public:
    using wavefrontFunctionType = double(*)(doubleThreeVector);
    using delWavefrontFunctionType = doubleThreeVector(*)(doubleThreeVector);
    using rabifreqFunctionType  = double(*)(doubleThreeVector, __float128);

    AISLaserBeam(doubleThreeVector k, __float128 omega,
                 wavefrontFunctionType wavefrontFunction,
                 delWavefrontFunctionType delWavefrontFunction,
                 rabifreqFunctionType rabiFreqFunction);
    ~AISLaserBeam();

    // At the moment I suppose that I can somehow pass a function as arguments
    // of the constructor to get the phases etc. Need to check on WiFi
    double GetPhi(doubleThreeVector pos); 
    doubleThreeVector GetDelPhi(doubleThreeVector pos);
    double GetRabiFreq(doubleThreeVector pos, __float128 t0, __float128 t);
    
    doubleThreeVector GetK();
    void SetK(doubleThreeVector k);

    __float128 GetOmega();
    void SetOmega(__float128 omega);

private:
    doubleThreeVector k;
    __float128 omega;

    wavefrontFunctionType wavefrontFunction;
    delWavefrontFunctionType delWavefrontFunction;
    rabifreqFunctionType rabifreqFunction;
};


#endif