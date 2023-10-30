#ifndef AISUTILITIES_HH
#define AISUTILITIES_HH

#include <iostream>
#include <array>
#include <cmath>
#include <string>
#include <complex>
#include <sstream>
#include <quadmath.h>
#include <iomanip>

/* Custom variable types*/
using doubleThreeVector = std::array<double, 3>;
using double3x3Matrix = std::array<std::array<double, 3>, 3>;

using quadThreeVector = std::array<__float128, 3>;
using quad3x3Matrix = std::array<std::array<__float128, 3>, 3>;

using complexDouble = std::complex<double>;

// complex number routines
double complexAbs(complexDouble z);
double complexArg(complexDouble z);

// Printing routine
void printMatrix(double3x3Matrix A);
void printVector(doubleThreeVector v);
void printQuadVariable(__float128 value);

/* Basic linear algebra functions*/
// Quad precision
// quad3x3Matrix matrixMultiply(const quad3x3Matrix& A, const quad3x3Matrix& B);
// quadThreeVector matrixMultiply(const quad3x3Matrix& A, const quadThreeVector& B);
// quad3x3Matrix matrixScalarMultiply(const quad3x3Matrix& A, const __float128& B);
// quad3x3Matrix matrixAdd(const quad3x3Matrix& A, const quad3x3Matrix& B);
// __float128 dotProduct(const quadThreeVector& A, const quadThreeVector& B);

// Double precision - multiplication
double3x3Matrix   dotProduct(const double3x3Matrix& A, const double3x3Matrix& B);
doubleThreeVector dotProduct(const doubleThreeVector& A, const double3x3Matrix& B);
doubleThreeVector dotProduct(const double3x3Matrix& A, const doubleThreeVector& B);
double            dotProduct(const doubleThreeVector& A, const doubleThreeVector& B);
double3x3Matrix   scalarMultiply(const double3x3Matrix& A, const double& B);
doubleThreeVector scalarMultiply(const doubleThreeVector& A, const double& B);
// Double precision - addition
double3x3Matrix   matrixAdd(const double3x3Matrix& A, const double3x3Matrix& B);
doubleThreeVector matrixAdd(const doubleThreeVector& A, const doubleThreeVector& B);
// Double precision - transposition
double3x3Matrix transpose(const double3x3Matrix& A);

/* Data type conversion routines*/
double convertScalarToDouble(__float128 x);
doubleThreeVector convertThreeVectorToDouble(quadThreeVector v);
double3x3Matrix convert3x3MatrixToDouble(quad3x3Matrix A);

#endif