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

    // Phase components which require quad precision
    __float128 GetPhaseQuad();
    void SetPhaseQuad(__float128 phase);
    // Phase components which do not require quad precision
    double GetPhaseDouble();
    void SetPhaseDouble(double phase);
    // Setter and getter for the phase error (double)
    double GetPhaseDoubleError();
    void SetPhaseDoubleError(double err);

    double GetAmplitude();
    void SetAmplitude(double amplitude);

    int GetState();
    void SetState(int state);

    int GetWavePacketID();
    void SetWavePacketID(int ID);

private:
    __float128 fTime = 0.0q;
    double fTime64 = 0.0;
    doubleThreeVector fPosition = {0.,0.,0.};
    doubleThreeVector fVelocity = {0.,0.,0.};
    __float128 fPhaseQuad = 0.0q;
    double fPhaseDouble = 0.0;
    double fPhaseDoubleError = 0.0;
    double fAmplitude = 1.0;
    int fState = 0;
    int fWavePacketID = 0;
};

#endif