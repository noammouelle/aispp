#ifndef AISWAVEPACKET_HH
#define AISWAVEPACKET_HH

#include "AISUtilities.hh"

class AISWavePacket{
public:
    AISWavePacket();
    ~AISWavePacket();

    __float128 GetTime();
    double GetTime64();
    void SetTime(__float128 time);

    doubleThreeVector GetPosition();
    void SetPosition(doubleThreeVector position);

    doubleThreeVector GetVelocity();
    void SetVelocity(doubleThreeVector velocity);

    __float128 GetPhase();
    void SetPhase(__float128 phase);

    __float128 GetAmplitude();
    void SetAmplitude(__float128 amplitude);

    int GetState();
    void SetState(int state);

    int GetWavePacketID();
    void SetWavePacketID(int ID);

private:
    __float128 fTime = 0.0q;
    double fTime64 = 0.0;
    doubleThreeVector fPosition = {0.,0.,0.};
    doubleThreeVector fVelocity = {0.,0.,0.};
    __float128 fPhase = 0.0q;
    __float128 fAmplitude = 1.0q;
    int fState = 0;
    int fWavePacketID = 0;
};

#endif