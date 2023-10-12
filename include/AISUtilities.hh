#ifndef AISUTILITIES_HH
#define AISUTILITIES_HH

#include <iostream>
#include <array>

/* Custom variable types*/
using doubleThreeVector = std::array<double, 3>;
using double3x3Matrix = std::array<std::array<double, 3>, 3>;

using quadThreeVector = std::array<__float128, 3>;
using quad3x3Matrix = std::array<std::array<__float128, 3>, 3>;

// Printing routine
void printMatrix(quad3x3Matrix A);
void printVector(quadThreeVector v);

/* Basic linear algebra functions*/
quad3x3Matrix matrixMultiply(const quad3x3Matrix& A, const quad3x3Matrix& B);
quadThreeVector matrixMultiply(const quad3x3Matrix& A, const quadThreeVector& B);
double3x3Matrix matrixMultiply(const double3x3Matrix& A, const double3x3Matrix& B);
doubleThreeVector matrixMultiply(const double3x3Matrix& A, const doubleThreeVector& B);
quad3x3Matrix matrixScalarMultiply(const quad3x3Matrix& A, const __float128& B);
__float128 dotProduct(const quadThreeVector& A, const quadThreeVector& B);
quad3x3Matrix matrixAdd(const quad3x3Matrix& A, const quad3x3Matrix& B);
double3x3Matrix matrixAdd(const double3x3Matrix& A, const double3x3Matrix& B);
doubleThreeVector matrixAdd(const doubleThreeVector& A, const doubleThreeVector& B);

/* Data type conversion routines*/
double convertScalarToDouble(__float128 x);
doubleThreeVector convertThreeVectorToDouble(quadThreeVector v);
double3x3Matrix convert3x3MatrixToDouble(quad3x3Matrix A);

#endif