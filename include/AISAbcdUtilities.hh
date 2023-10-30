#ifndef AISABCDUTILITIES_HH
#define AISABCDUTILITIES_HH

#include "AISUtilities.hh"

// define functions to compute the abcd matrices
/* 
    In previous implementations, ABCD are expressed as a function
    of t0 and t1, i.e.
        A(t0, t1)
    but really they are functions of the time interval T = t1 - t0. 
    This is how T should be defined here.
*/
double3x3Matrix computeA(const double3x3Matrix& alpha, const double3x3Matrix& gamma, 
                             const doubleThreeVector& gVector, const double& T);

double3x3Matrix computeB(const double3x3Matrix& alpha, const double3x3Matrix& gamma, 
                             const doubleThreeVector& gVector, const double& T);

double3x3Matrix computeC(const double3x3Matrix& alpha, const double3x3Matrix& gamma, 
                             const doubleThreeVector& gVector, const double& T);

double3x3Matrix computeD(const double3x3Matrix& alpha, const double3x3Matrix& gamma, 
                             const doubleThreeVector& gVector, const double& T);

doubleThreeVector computeXi(const double3x3Matrix& alpha, const double3x3Matrix& gamma, 
                             const doubleThreeVector& gVector, const double& T);

doubleThreeVector computePhi(const double3x3Matrix& alpha, const double3x3Matrix& gamma, 
                             const doubleThreeVector& gVector, const double& T);

#endif