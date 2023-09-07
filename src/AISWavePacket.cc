#include "AISWavePacket.hh"

AISWavePacket::AISWavePacket(){};
AISWavePacket::~AISWavePacket(){};

__float128 AISWavePacket::GetTime(){
    return fTime;
};
double AISWavePacket::GetTime64(){
    return fTime64;
};
void AISWavePacket::SetTime(__float128 time){
    fTime = time;
    fTime64 = time;
};

doubleThreeVector AISWavePacket::GetPosition(){
    return fPosition;
};
void AISWavePacket::SetPosition(doubleThreeVector position){
    fPosition = position;
};

doubleThreeVector AISWavePacket::GetVelocity(){
    return fVelocity;
};
void AISWavePacket::SetVelocity(doubleThreeVector velocity){
    fVelocity = velocity;
};

__float128 AISWavePacket::GetPhase(){
    return fPhase;
};
void AISWavePacket::SetPhase(__float128 phase){
    fPhase = phase;
}

__float128 AISWavePacket::GetAmplitude(){
    return fAmplitude;
};
void AISWavePacket::SetAmplitude(__float128 amplitude){
    fAmplitude = amplitude;
};

int AISWavePacket::GetState(){
    return fState;
};
void AISWavePacket::SetState(int state){
    fState = state;
};

int AISWavePacket::GetWavePacketID(){
    return fWavePacketID;
};
void AISWavePacket::SetWavePacketID(int ID){
    fWavePacketID = ID;
}