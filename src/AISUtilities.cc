#include "AISUtilities.hh"
#include <iostream>

// complex number routines
double complexAbs(complexDouble z){
    double result = std::abs(z);
    return result;
}

double complexArg(complexDouble z){
    double result = std::arg(z);
    return result;
}

// Printing routine
void printMatrix(double3x3Matrix A){
    for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                std::cout << A[i][j] << "\t"; 
            }
            std::cout << std::endl; 
        }
    std::cout << "" <<std::endl;
}

void printVector(doubleThreeVector v){
    for(int i = 0; i < 3; ++i){
        std::cout << v[i] << std::endl;
    }
    std::cout << "" <<std::endl;
}

void printQuadVariable(__float128 value) {
    // Convert the __float128 to long double
    long double ldValue = static_cast<long double>(value);

    // Print the long double with full precision
    std::cout << std::setprecision(34) << ldValue << std::endl;
}

/* Basic linear algebra functions*/
// Quad precision
quad3x3Matrix matrixMultiply(const quad3x3Matrix& A, const quad3x3Matrix& B){

    quad3x3Matrix result = {};

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}

quadThreeVector matrixMultiply(const quad3x3Matrix& A, const quadThreeVector& B){

    quadThreeVector result = {};

    for (int i = 0; i < 3; ++i) {
        for (int k = 0; k < 3; ++k) {
            result[i] += A[i][k] * B[k];
        }
    }
    return result;
}

quad3x3Matrix matrixScalarMultiply(const quad3x3Matrix& A, const __float128& B){

    quad3x3Matrix result = {};

    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < 3; ++j){
            result[i][j] *= B;
        }
    }
    return result;
}

__float128 dotProduct(const quadThreeVector& A, const quadThreeVector& B){

    double result = 0.0;

    for (int i = 0; i < 3; ++i) {
        result += A[i] * B[i];
    }
    return result;
}

quad3x3Matrix matrixAdd(const quad3x3Matrix& A, const quad3x3Matrix& B){
    
    quad3x3Matrix result = {};

    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < 3; ++j){
            result[i][j] = A[i][j] + B[i][j];
        }
    }
    return result;
}

// Double precision - multiplication
double3x3Matrix dotProduct(const double3x3Matrix& A, const double3x3Matrix& B){

    double3x3Matrix result = {};

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}

doubleThreeVector dotProduct(const double3x3Matrix& A, const doubleThreeVector& B){

    doubleThreeVector result = {};

    for (int i = 0; i < 3; ++i) {
        for (int k = 0; k < 3; ++k) {
            result[i] += A[i][k] * B[k];
        }
    }
    return result;
}

doubleThreeVector dotProduct(const doubleThreeVector& A, const double3x3Matrix& B){

    doubleThreeVector result = {};

    for (int i = 0; i < 3; ++i) {
        for (int k = 0; k < 3; ++k) {
            result[i] += A[k] * B[k][i];
        }
    }
    return result;
}

double dotProduct(const doubleThreeVector& A, const doubleThreeVector& B){

    double result = 0.0;

    for (int i = 0; i < 3; ++i) {
        result += A[i] * B[i];
    }
    return result;
}

double3x3Matrix scalarMultiply(const double3x3Matrix& A, const double& B){
    
    double3x3Matrix result = {};

    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < 3; ++j){ 
            result[i][j] = A[i][j] * B;
        }
    }
    return result;
}

doubleThreeVector scalarMultiply(const doubleThreeVector& A, const double& B){
    
    doubleThreeVector result = {};

    for(int i = 0; i < 3; ++i){
        result[i] = A[i] * B;
    }
    return result;
}

double3x3Matrix matrixAdd(const double3x3Matrix& A, const double3x3Matrix& B){
    
    double3x3Matrix result = {};

    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < 3; ++j){
            result[i][j] = A[i][j] + B[i][j];
        }
    }
    return result;
}

doubleThreeVector matrixAdd(const doubleThreeVector& A, const doubleThreeVector& B){
    
    doubleThreeVector result = {};

    for(int i = 0; i < 3; ++i){
        result[i] = A[i] + B[i];
    }
    return result;
}

// Double precision - transposition
double3x3Matrix transpose(const double3x3Matrix& A){
    
    double3x3Matrix result = {};

    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < 3; ++j){ 
            result[i][j] = A[j][i];
        }
    }
    return result;
}

/* Data type conversion routines*/
double convertScalarToDouble(__float128 x){
    double result = static_cast<double>(x);
    return result;
}

doubleThreeVector convertThreeVectorToDouble(quadThreeVector v){
    doubleThreeVector result;
    for(int i = 0; i < 3; ++i){
        result[i] = static_cast<double>(v[i]);
    }
    return result;
}

double3x3Matrix convert3x3MatrixToDouble(quad3x3Matrix A){
    double3x3Matrix result;
    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < 3; ++j){
            result[i][j] = static_cast<double>(A[i][j]);
        }
    }
    return result;
}

doubleThreeVector computeDetunedWaveVector(const doubleThreeVector& k0, const double& vz)
{
    doubleThreeVector detunedK = k0;
    for(int i = 0; i < 3; ++i)
    {
        if(i == 2){
            detunedK[i] = k0[i] /(1 - vz / c);
        }
    }

    return detunedK;
}

__float128 computeDetunedOmega(const __float128& omega0, const double& vz){
    return omega0 / (1 - vz / c);
}

void writeAtomEnsembleToFile(std::string fName, AISAtomEnsemble* atomEnsemble)
{
    // write a file to cross check
    std::ofstream outFile(fName);
    std::cout << "Writing to file: " << fName << std::endl;
    outFile << std::fixed; // Use fixed-point notation
    outFile << std::setprecision(15);

    outFile << "State, Amplitude, X, Y, Z, VX, VY, VZ, Phase, PhaseQuad" << std::endl;
    for(int i = 0; i < atomEnsemble->GetNumberOfAtoms(); ++i)
    {
        AISAtom* currentAtom   = atomEnsemble->GetAtom(i);
        for(int j = 0; j < currentAtom->GetNumberOfWavePackets(); ++j)
        {
            AISWavePacket* currentWavePacket = currentAtom->GetWavePacket(j);
            int currentState = currentWavePacket->GetState();
            double currentAmplitude = currentWavePacket->GetAmplitude();
            doubleThreeVector currentPos = currentWavePacket->GetPosition();
            doubleThreeVector currentVel = currentWavePacket->GetVelocity();
            double currentPhaseDouble = currentWavePacket->GetPhaseDouble();
            __float128 currentPhaseQuad = currentWavePacket->GetPhaseQuad();

            long double currentPhaseQuadAsLongDouble = static_cast<long double>(currentPhaseQuad);

            char currentPhaseQuadStr[50];  // Adjust the buffer size as needed
            snprintf(currentPhaseQuadStr, sizeof(currentPhaseQuadStr), "%.34Le", currentPhaseQuadAsLongDouble);

            outFile << currentState << "," << currentAmplitude << "," << currentPos[0] << "," << currentPos[1] << "," << currentPos[2] << "," 
                    << currentVel[0] << "," << currentVel[1] << "," << currentVel[2] << ","
                    << currentPhaseDouble << "," << currentPhaseQuadStr << std::endl;
        }
    }

    outFile.close();
}