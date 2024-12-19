#ifndef AISLASERBEAM_HH
#define AISLASERBEAM_HH

#include "AISConstants.hh"
#include "AISUtilities.hh"

#include <memory>

class AISLaserBeam
{
public:
    using wavefrontFunctionType = double(*)(const doubleThreeVector&);
    using delWavefrontFunctionType = doubleThreeVector(*)(const doubleThreeVector&);
    using rabifreqFunctionType  = double(*)(const doubleThreeVector&, const __float128&, const __float128&);

    AISLaserBeam(doubleThreeVector k, __float128 omega, double rabiFreq, double phi0,
                 std::shared_ptr<wavefrontFunctionType> wavefrontFunction,
                 std::shared_ptr<delWavefrontFunctionType> delWavefrontFunction,
                 std::shared_ptr<rabifreqFunctionType> rabiFreqFunction);
    ~AISLaserBeam();

    // At the moment I suppose that I can somehow pass a function as arguments
    // of the constructor to get the phases etc. Need to check on WiFi
    double GetPhi(const doubleThreeVector& pos); 
    doubleThreeVector GetDelPhi(const doubleThreeVector& pos);
    double GetRabiFreq(const doubleThreeVector& pos, const __float128& t0, const __float128& t);
    
    doubleThreeVector GetK(__float128 t);
    void SetK(doubleThreeVector k);

    __float128 GetOmega(__float128 t);
    void SetOmega(__float128 omega);

    __float128 GetFrequencyChirp();
    void SetFrequencyChirp(__float128 frequencyChirp);

    doubleThreeVector GetKChirp();
    void SetKChirp(doubleThreeVector kChirp);

private:
    double rabiFreq; // central rabi frequency

    doubleThreeVector k;
    __float128 omega, frequencyChirp;
    doubleThreeVector kChirp;
    double phi0;

    std::shared_ptr<wavefrontFunctionType> wavefrontFunction;
    std::shared_ptr<delWavefrontFunctionType> delWavefrontFunction;
    std::shared_ptr<rabifreqFunctionType> rabifreqFunction;
};


#endif