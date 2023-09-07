#ifndef AISWAVEFRONT_HH
#define AISWAVEFRONT_HH

#include <array>

using threeVector = std::array<double, 3>;

class AISWaveFront
{
private:
    /* data */
public:
    AISWaveFront(/* args */);
    ~AISWaveFront();

    void GetValue(threeVector pos);
};

#endif