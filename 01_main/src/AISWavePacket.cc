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
    fTime64 = convertScalarToDouble(time);
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

__float128 AISWavePacket::GetPhaseQuad(){
    return fPhaseQuad;
};
void AISWavePacket::SetPhaseQuad(__float128 phase){
    fPhaseQuad = phase;
}

double AISWavePacket::GetPhaseDouble(){
    return fPhaseDouble;
};
void AISWavePacket::SetPhaseDouble(double phase){
    fPhaseDouble = phase;
}

double AISWavePacket::GetPhaseDoubleError(){
    return fPhaseDoubleError;
};
void AISWavePacket::SetPhaseDoubleError(double err){
    fPhaseDoubleError = err;
}

double AISWavePacket::GetAmplitude(){
    return fAmplitude;
};
void AISWavePacket::SetAmplitude(double amplitude){
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

std::string AISWavePacket::GetPath(){
    return path;
};
void AISWavePacket::SetPath(std::string aPath){
    path = aPath;
};

bool AISWavePacket::GetWillInterfere(){
    return willInterfere;
};
void AISWavePacket::SetWillInterfere(bool willInterfereVal){
    willInterfere = willInterfereVal;
};

doubleThreeVector AISWavePacket::GetPosStar(){
    return posStar;
};

void AISWavePacket::SetPosStar(doubleThreeVector posStarVal){
    posStar = posStarVal;
};

doubleThreeVector AISWavePacket::GetVelStar(){
    return velStar;
};

void AISWavePacket::SetVelStar(doubleThreeVector velStarVal){
    velStar = velStarVal;
};

std::vector<std::string> AISWavePacket::GetDetectablePaths(){
    return detectablePaths;
};

void AISWavePacket::SetDetectablePaths(std::vector<std::string> detectablePathsVal){
    detectablePaths = detectablePathsVal;
};

doubleThreeVector AISWavePacket::GetPos0(){
    return pos0;
};

void AISWavePacket::SetPos0(doubleThreeVector pos0Val){
    pos0 = pos0Val;
};

doubleThreeVector AISWavePacket::GetVel0(){
    return vel0;
};

void AISWavePacket::SetVel0(doubleThreeVector vel0Val){
    vel0 = vel0Val;
};

__float128 AISWavePacket::GetT0(){
    return t0;
};

void AISWavePacket::SetT0(__float128 t0Val){
    t0 = t0Val;
};