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
// all in .hh now for inlining.


// Double precision - addition
// all in .hh now for inlining.

// Double precision - transposition
// all in .hh now for inlining.

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

// Function to convert std::array to GSL matrix
void arrayToGslMatrix(const double6x6Matrix& src, gsl_matrix* dst) {
    for (size_t i = 0; i < 6; ++i) {
        for (size_t j = 0; j < 6; ++j) {
            gsl_matrix_set(dst, i, j, src[i][j]);
        }
    }
}

// Function to convert GSL matrix back to std::array
void gslMatrixToArray(const gsl_matrix* src, double6x6Matrix& dst) {
    for (size_t i = 0; i < 6; ++i) {
        for (size_t j = 0; j < 6; ++j) {
            dst[i][j] = gsl_matrix_get(src, i, j);
        }
    }
}

// Function to compute the inverse of a 6x6 matrix
double6x6Matrix invertMatrix(const double6x6Matrix& M) {
    double6x6Matrix inverseMatrix;

    // Create GSL matrices for the original and inverse
    gsl_matrix* gslM = gsl_matrix_alloc(6, 6);
    gsl_matrix* gslInverse = gsl_matrix_alloc(6, 6);
    gsl_permutation* perm = gsl_permutation_alloc(6);
    int signum;

    // Copy data from M to the GSL matrix
    arrayToGslMatrix(M, gslM);

    // Perform LU decomposition
    gsl_linalg_LU_decomp(gslM, perm, &signum);

    // Compute the inverse
    gsl_linalg_LU_invert(gslM, perm, gslInverse);

    // Copy the result back to the std::array
    gslMatrixToArray(gslInverse, inverseMatrix);

    // Free GSL resources
    gsl_matrix_free(gslM);
    gsl_matrix_free(gslInverse);
    gsl_permutation_free(perm);

    return inverseMatrix;
}

// Function to matrix exponential of a 6x6 matrix
double6x6Matrix expMatrix(const double6x6Matrix& M) {
    double6x6Matrix expM;

    // Create GSL matrices for the original and inverse
    gsl_matrix* gslM = gsl_matrix_alloc(6, 6);
    gsl_matrix* gslExpM = gsl_matrix_alloc(6, 6);

    // Copy data from M to the GSL matrix
    arrayToGslMatrix(M, gslM);

    // Compute the exponential
    gsl_linalg_exponential_ss(gslM, gslExpM, GSL_PREC_DOUBLE);

    // Copy the result back to the std::array
    gslMatrixToArray(gslExpM, expM);

    // Free GSL resources
    gsl_matrix_free(gslM);
    gsl_matrix_free(gslExpM);

    return expM;
}

std::string float128ToString(__float128 value) {
    char buffer[128];
    quadmath_snprintf(buffer, sizeof(buffer), "%.36Qg", value);
    return std::string(buffer);
}

// Function to compute the factorial of an integer
int factorial(int n) {
    if (n == 0) {
        return 1;
    }
    return n * factorial(n - 1);
}

// Function to convert Noll index to Zernike indices
std::array<int, 2> nollToZernike(int j) {
    // from https://github.com/opticspy/lightpipes/blob/master/LightPipes/zernike.py
    if (j == 0) {
        // Noll indices start at 1 so 0 is invalid
        throw std::invalid_argument("Noll index must be greater than 0");
    }
    
    int n = 0;
    int j1 = j-1;
    while (j1 > n) {
        n++;
        j1 -= n;
    }
    int m = pow(-1,j) * ((n % 2) + 2 * std::floor((j1+((n+1)%2)) / 2.0 ));

    return {n, m};
}