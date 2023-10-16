#include "AISUtilities.hh"
#include <iostream>

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