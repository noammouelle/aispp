#ifndef AISUTILITIES_HH
#define AISUTILITIES_HH

#include "AISConstants.hh"

#include <iostream>
#include <array>
#include <vector>
#include <cmath>
#include <string>
#include <complex>
#include <sstream>
#include <quadmath.h>
#include <iomanip>
#include <iomanip>
#include <fstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

/* Custom variable types*/
using doubleThreeVector = std::array<double, 3>;
using double3x3Matrix = std::array<std::array<double, 3>, 3>;

using doubleSixVector = std::array<double,6>;
using double6x6Matrix = std::array<std::array<double,6>,6>;

using quadThreeVector = std::array<__float128, 3>;
using quad3x3Matrix = std::array<std::array<__float128, 3>, 3>;

using threeVectorList = std::vector<doubleThreeVector>;
using sixVectorList = std::vector<std::array<double, 6>>;

using complexDouble = std::complex<double>;

using intVector = std::vector<int>;
using intTuple = std::array<intVector,2>;
using doubleVector = std::vector<double>;
using quadVector = std::vector<__float128>;
using boolVector = std::vector<bool>;

using wavefrontFunctionType = double(*)(const doubleThreeVector&, const double&, const double&);
using delWavefrontFunctionType = doubleThreeVector(*)(const doubleThreeVector&, const double&, const double&);
using rabifreqFunctionType  = double(*)(const doubleThreeVector&, const double&, const __float128&,const  __float128&);

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
doubleSixVector   dotProduct(const double6x6Matrix& A, const doubleSixVector& b);
doubleThreeVector dotProduct(const doubleThreeVector& A, const double3x3Matrix& B);
doubleThreeVector dotProduct(const double3x3Matrix& A, const doubleThreeVector& B);
double            dotProduct(const doubleThreeVector& A, const doubleThreeVector& B);
double3x3Matrix   scalarMultiply(const double3x3Matrix& A, const double& B);
double6x6Matrix   scalarMultiply(const double6x6Matrix& A, const double& B);
doubleThreeVector scalarMultiply(const doubleThreeVector& A, const double& B);
// Double precision - addition
double3x3Matrix   matrixAdd(const double3x3Matrix& A, const double3x3Matrix& B);
double6x6Matrix   matrixAdd(const double6x6Matrix& A, const double6x6Matrix& B);
doubleSixVector   matrixAdd(const doubleSixVector& A, const doubleSixVector& B);
doubleThreeVector matrixAdd(const doubleThreeVector& A, const doubleThreeVector& B);
// Double precision - transposition
double3x3Matrix transpose(const double3x3Matrix& A);

/* Data type conversion routines*/
double convertScalarToDouble(__float128 x);
doubleThreeVector convertThreeVectorToDouble(quadThreeVector v);
double3x3Matrix convert3x3MatrixToDouble(quad3x3Matrix A);

/* Detuning functions */
doubleThreeVector computeDetunedWaveVector(const doubleThreeVector& k0, const double& vz);
__float128 computeDetunedOmega(const __float128& omega0, const double& vz);

/*Linalg and conversion routines*/
void arrayToGslMatrix(const double6x6Matrix& src, gsl_matrix* dst);
void gslMatrixToArray(const gsl_matrix* src, double6x6Matrix& dst);
double6x6Matrix invertMatrix(const double6x6Matrix& M);
double6x6Matrix expMatrix(const double6x6Matrix& M);

/*Printing routines*/
std::string float128ToString(__float128 value);

/*Math routines*/
// factorial of an integer
int factorial(int n);

/*Zernike routines*/
std::array<int,2> nollToZernike(int noll);



#endif