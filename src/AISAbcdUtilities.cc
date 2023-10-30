#include "AISAbcdUtilities.hh"

double3x3Matrix computeA(const double3x3Matrix& alpha, const double3x3Matrix& gamma, 
                             const doubleThreeVector& gVector, const double& T)
{
    double3x3Matrix A;
    double3x3Matrix I = {{ {{1.0, 0.0, 0.0}},
                           {{0.0, 1.0, 0.0}},
                           {{0.0, 0.0, 1.0}} }};

    // Precompute all the required coefficients
    double3x3Matrix gamma2 = dotProduct(gamma, gamma);

    double3x3Matrix alpha2 = dotProduct(alpha, alpha);
    double3x3Matrix alpha3 = dotProduct(alpha2, alpha);
    double3x3Matrix alpha4 = dotProduct(alpha3,alpha);

    double3x3Matrix alphaDotGamma = dotProduct(alpha, gamma);
    double3x3Matrix gammaDotAlpha = dotProduct(gamma, alpha);

    double3x3Matrix alpha2DotGamma = dotProduct(alpha2, gamma);
    double3x3Matrix gammaDotAlpha2 = dotProduct(gamma, alpha2);

    double3x3Matrix alphaDotGammaDotAlpha = dotProduct(alphaDotGamma, alpha);

    double T2 = T  * T;
    double T3 = T2 * T;
    double T4 = T3 * T;
    double T5 = T4 * T;
    double T6 = T5 * T;

    /* Compute the Taylor series */
    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < 3; ++j){
            A[i][j] = I[i][j] + alpha[i][j] * T
                    + (alpha2[i][j] + gamma[i][j])*T2 / 2
                    + (alpha3[i][j] + 2*alphaDotGamma[i][j] + gammaDotAlpha[i][j])*T3 / 6
                    + (alpha4[i][j] + 3*alpha2DotGamma[i][j] + 2*alphaDotGammaDotAlpha[i][j]
                        + gammaDotAlpha2[i][j] + gamma2[i][j]) * T4 / 24 ;
        }
    }

    return A;
}

double3x3Matrix computeB(const double3x3Matrix& alpha, const double3x3Matrix& gamma, 
                             const doubleThreeVector& gVector, const double& T)
{
    double3x3Matrix B;

    double3x3Matrix I = {{ {{1.0, 0.0, 0.0}},
                           {{0.0, 1.0, 0.0}},
                           {{0.0, 0.0, 1.0}} }};

    // Precompute all the required coefficients
    double3x3Matrix gamma2 = dotProduct(gamma, gamma);

    double3x3Matrix alpha2 = dotProduct(alpha, alpha);
    double3x3Matrix alpha3 = dotProduct(alpha2, alpha);
    double3x3Matrix alpha4 = dotProduct(alpha3,alpha);

    double3x3Matrix alphaDotGamma = dotProduct(alpha, gamma);
    double3x3Matrix gammaDotAlpha = dotProduct(gamma, alpha);

    double3x3Matrix alpha2DotGamma = dotProduct(alpha2, gamma);
    double3x3Matrix gammaDotAlpha2 = dotProduct(gamma, alpha2);

    double3x3Matrix alphaDotGammaDotAlpha = dotProduct(alphaDotGamma, alpha);

    double T2 = T  * T;
    double T3 = T2 * T;
    double T4 = T3 * T;
    double T5 = T4 * T;
    double T6 = T5 * T;

    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < 3; ++j){
            B[i][j] = I[i][j] * T + alpha[i][j] * T2 
                    + (3 * alpha2[i][j] + gamma[i][j])*T3 / 6
                    + (5 * alpha4[i][j] + 3*alpha2DotGamma[i][j] + 4*alphaDotGammaDotAlpha[i][j]
                        + 3 * gammaDotAlpha2[i][j] + gamma2[i][j]) * T5 / 120 ;
        }
    }

    return B;
}

double3x3Matrix computeC(const double3x3Matrix& alpha, const double3x3Matrix& gamma, 
                             const doubleThreeVector& gVector, const double& T)
{
    double3x3Matrix C;

    double3x3Matrix I = {{ {{1.0, 0.0, 0.0}},
                           {{0.0, 1.0, 0.0}},
                           {{0.0, 0.0, 1.0}} }};

    // Precompute all the required coefficients
    double3x3Matrix gamma2 = dotProduct(gamma, gamma);

    double3x3Matrix alpha2 = dotProduct(alpha, alpha);
    double3x3Matrix alpha3 = dotProduct(alpha2, alpha);
    double3x3Matrix alpha4 = dotProduct(alpha3,alpha);

    double3x3Matrix alphaDotGamma = dotProduct(alpha, gamma);
    double3x3Matrix gammaDotAlpha = dotProduct(gamma, alpha);

    double3x3Matrix alpha2DotGamma = dotProduct(alpha2, gamma);
    double3x3Matrix gammaDotAlpha2 = dotProduct(gamma, alpha2);

    double3x3Matrix alphaDotGammaDotAlpha = dotProduct(alphaDotGamma, alpha);

    double T2 = T  * T;
    double T3 = T2 * T;
    double T4 = T3 * T;
    double T5 = T4 * T;
    double T6 = T5 * T;

    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < 3; ++j){
            C[i][j] = gamma[i][j] * T
                    + (alphaDotGamma[i][j] + gammaDotAlpha[i][j]) * T2 / 2
                    + (alpha2DotGamma[i][j] + alphaDotGammaDotAlpha[i][j]
                        + gammaDotAlpha2[i][j] + gamma2[i][j]) * T3 / 6 ;
        }
    }

    return C;
}

double3x3Matrix computeD(const double3x3Matrix& alpha, const double3x3Matrix& gamma, 
                             const doubleThreeVector& gVector, const double& T)
{
    double3x3Matrix D;

    double3x3Matrix I = {{ {{1.0, 0.0, 0.0}},
                           {{0.0, 1.0, 0.0}},
                           {{0.0, 0.0, 1.0}} }};

    // Precompute all the required coefficients
    double3x3Matrix gamma2 = dotProduct(gamma, gamma);

    double3x3Matrix alpha2 = dotProduct(alpha, alpha);
    double3x3Matrix alpha3 = dotProduct(alpha2, alpha);
    double3x3Matrix alpha4 = dotProduct(alpha3,alpha);

    double3x3Matrix alphaDotGamma = dotProduct(alpha, gamma);
    double3x3Matrix gammaDotAlpha = dotProduct(gamma, alpha);

    double3x3Matrix alpha2DotGamma = dotProduct(alpha2, gamma);
    double3x3Matrix gammaDotAlpha2 = dotProduct(gamma, alpha2);

    double3x3Matrix alphaDotGammaDotAlpha = dotProduct(alphaDotGamma, alpha);

    double T2 = T  * T;
    double T3 = T2 * T;
    double T4 = T3 * T;
    double T5 = T4 * T;
    double T6 = T5 * T;

    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < 3; ++j){
            D[i][j] = I[i][j] + alpha[i][j] * T 
                    + (alpha2[i][j] + gamma[i][j]) * T2 / 2
                    + (alpha3[i][j] + alphaDotGamma[i][j] + 2 * gammaDotAlpha[i][j]) * T3 / 6
                    + (alpha4[i][j] + alpha2DotGamma[i][j]
                        + 2 * alphaDotGammaDotAlpha[i][j]
                        + 3 * gammaDotAlpha2[i][j]
                        + gamma2[i][j]) * T4 / 24 ;
        }
    }

    return D;
}

doubleThreeVector computeXi(const double3x3Matrix& alpha, const double3x3Matrix& gamma, 
                             const doubleThreeVector& gVector, const double& T)
{
    doubleThreeVector Xi;
    double3x3Matrix xiMatrix;

    double3x3Matrix I = {{ {{1.0, 0.0, 0.0}},
                           {{0.0, 1.0, 0.0}},
                           {{0.0, 0.0, 1.0}} }};

    // Precompute all the required coefficients
    double3x3Matrix gamma2 = dotProduct(gamma, gamma);

    double3x3Matrix alpha2 = dotProduct(alpha, alpha);
    double3x3Matrix alpha3 = dotProduct(alpha2, alpha);
    double3x3Matrix alpha4 = dotProduct(alpha3,alpha);

    double3x3Matrix alphaDotGamma = dotProduct(alpha, gamma);
    double3x3Matrix gammaDotAlpha = dotProduct(gamma, alpha);

    double3x3Matrix alpha2DotGamma = dotProduct(alpha2, gamma);
    double3x3Matrix gammaDotAlpha2 = dotProduct(gamma, alpha2);

    double3x3Matrix alphaDotGammaDotAlpha = dotProduct(alphaDotGamma, alpha);

    double T2 = T  * T;
    double T3 = T2 * T;
    double T4 = T3 * T;
    double T5 = T4 * T;
    double T6 = T5 * T;

    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < 3; ++j){
            xiMatrix[i][j] = I[i][j] * T2 / 2 + alpha[i][j] * T3 / 6 
                            + (3 * alpha2[i][j] + gamma[i][j]) * T4 / 24
                            + (4 * alpha3[i][j] + 2 * alphaDotGamma[i][j] + 2 * gammaDotAlpha[i][j]) * T5 / 120
                            + (5 * alpha4[i][j] + 3 * alpha2DotGamma[i][j]
                                + 4 * alphaDotGammaDotAlpha[i][j]
                                + 3 * gammaDotAlpha2[i][j]
                                + gamma2[i][j]) * T6 / 720 ;
        }
    }
    Xi = dotProduct(xiMatrix, gVector);
    return Xi;
}

doubleThreeVector computePhi(const double3x3Matrix& alpha, const double3x3Matrix& gamma, 
                             const doubleThreeVector& gVector, const double& T)
{
    doubleThreeVector Phi;
    double3x3Matrix phiMatrix;

    double3x3Matrix I = {{ {{1.0, 0.0, 0.0}},
                           {{0.0, 1.0, 0.0}},
                           {{0.0, 0.0, 1.0}} }};

    // Precompute all the required coefficients
    double3x3Matrix gamma2 = dotProduct(gamma, gamma);

    double3x3Matrix alpha2 = dotProduct(alpha, alpha);
    double3x3Matrix alpha3 = dotProduct(alpha2, alpha);
    double3x3Matrix alpha4 = dotProduct(alpha3,alpha);

    double3x3Matrix alphaDotGamma = dotProduct(alpha, gamma);
    double3x3Matrix gammaDotAlpha = dotProduct(gamma, alpha);

    double3x3Matrix alpha2DotGamma = dotProduct(alpha2, gamma);
    double3x3Matrix gammaDotAlpha2 = dotProduct(gamma, alpha2);

    double3x3Matrix alphaDotGammaDotAlpha = dotProduct(alphaDotGamma, alpha);

    double T2 = T  * T;
    double T3 = T2 * T;
    double T4 = T3 * T;
    double T5 = T4 * T;
    double T6 = T5 * T;

    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < 3; ++j){
            phiMatrix[i][j] = I[i][j] * T + alpha[i][j] * T2 / 2 
                            + (alpha2[i][j] + gamma[i][j]) * T3 / 6        
                            + (alpha3[i][j] + alphaDotGamma[i][j] + 2 * gammaDotAlpha[i][j]) * T4 / 24
                            + (alpha4[i][j] + alpha2DotGamma[i][j] + 2 * alphaDotGammaDotAlpha[i][j]
                                + 3 * gammaDotAlpha2[i][j] + gamma2[i][j]) * T5 / 120 ; 
        }
    }
    Phi = dotProduct(phiMatrix, gVector);
    return Phi;
}