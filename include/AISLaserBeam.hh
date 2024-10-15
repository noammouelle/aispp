#ifndef AISLASERBEAM_HH
#define AISLASERBEAM_HH

#include "AISConstants.hh"
#include "AISUtilities.hh"

class AISLaserBeam
{
private:
    doubleThreeVector k;
    __float128 omega;

    AISIntensityProfile* intensityProfile;
    AISWavefront* wavefront;

public:
    AISLaserBeam(doubleThreeVector k, __float128 omega, 
                 AISIntensityProfile* intensityProfile,
                 AISWavefront* wavefront);
    ~AISLaserBeam();

    double GetPhase(doubleThreeVector pos); 
    double GetRabiFreq(doubleThreeVector pos, __float128 t);
    
    doubleThreeVector GetK();
    void SetK(doubleThreeVector k);

    __float128 GetOmega();
    void SetOmega(__float128 omega);
};


#endif