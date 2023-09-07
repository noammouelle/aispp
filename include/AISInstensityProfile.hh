#ifndef AISINTENSITYPROFILE_HH
#define AISINTENSITYPROFILE_HH

#include <array>

using threeVector = std::array<double,3>;

class AISInstensityProfile
{
private:
    /* data */
public:
    AISInstensityProfile(/* args */);
    ~AISInstensityProfile();

    double GetEffectiveRabiFreq(threeVector pos);
};



#endif