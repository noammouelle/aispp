#include "AISConstants.hh"
#include "AISAbcdPropagator.hh"


AISAbcdPropagator::AISAbcdPropagator(__float128 dt){
    deltaTime   = dt;
    deltaTime64 = convertScalarToDouble(deltaTime);

    // Identity
     I = {{ {{1.0q, 0.0q, 0.0q}},
            {{0.0q, 1.0q, 0.0q}},
            {{0.0q, 0.0q, 1.0q}} }};
    
    // Dummy values for coefficients
    alpha   = {};
    gamma   = {};
    gVector = {};
    
    // Compute the ABCD matrices
    setABCDXiPhi();
}

AISAbcdPropagator::~AISAbcdPropagator(){}

void AISAbcdPropagator::setABCDXiPhi(){

    // Precompute all the required coefficients
    quad3x3Matrix gamma2 = matrixMultiply(gamma, gamma);

    quad3x3Matrix alpha2 = matrixMultiply(alpha, alpha);
    quad3x3Matrix alpha3 = matrixMultiply(alpha2, alpha);
    quad3x3Matrix alpha4 = matrixMultiply(alpha3,alpha);

    quad3x3Matrix alphaDotGamma = matrixMultiply(alpha, gamma);
    quad3x3Matrix gammaDotAlpha = matrixMultiply(gamma, alpha);

    quad3x3Matrix alpha2DotGamma = matrixMultiply(alpha2, gamma);
    quad3x3Matrix gammaDotAlpha2 = matrixMultiply(gamma, alpha2);

    quad3x3Matrix alphaDotGammaDotAlpha = matrixMultiply(alphaDotGamma, alpha);

    __float128 dt  = deltaTime;
    __float128 dt2 = dt * dt;
    __float128 dt3 = dt2 * dt;
    __float128 dt4 = dt3 * dt;
    __float128 dt5 = dt4 * dt;
    __float128 dt6 = dt5 * dt;

    /* Compute the Taylor series */

    // A
    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < 3; ++j){
            A[i][j] = I[i][j] + alpha[i][j] * dt
                    + (alpha2[i][j] + gamma[i][j])*dt2 / 2
                    + (alpha3[i][j] + 2*alphaDotGamma[i][j] + gammaDotAlpha[i][j])*dt3 / 6
                    + (alpha4[i][j] + 3*alpha2DotGamma[i][j] + 2*alphaDotGammaDotAlpha[i][j]
                        + gammaDotAlpha2[i][j] + gamma2[i][j]) * dt4 / 24 ;
        }
    }

    // B
    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < 3; ++j){
            B[i][j] = I[i][j] * dt + alpha[i][j] * dt2 
                    + (3 * alpha2[i][j] + gamma[i][j])*dt3 / 6
                    + (5 * alpha4[i][j] + 3*alpha2DotGamma[i][j] + 4*alphaDotGammaDotAlpha[i][j]
                        + 3 * gammaDotAlpha2[i][j] + gamma2[i][j]) * dt5 / 120 ;
        }
    }

    // C
    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < 3; ++j){
            C[i][j] = gamma[i][j] * dt 
                    + (alphaDotGamma[i][j] + gammaDotAlpha[i][j]) * dt2 / 2
                    + (alpha2DotGamma[i][j] + alphaDotGammaDotAlpha[i][j]
                        + gammaDotAlpha2[i][j] + gamma2[i][j]) * dt3 / 6 ;
        }
    }

    // D
    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < 3; ++j){
            D[i][j] = I[i][j] + alpha[i][j] * dt 
                    + (alpha2[i][j] + gamma[i][j]) * dt2 / 2
                    + (alpha3[i][j] + alphaDotGamma[i][j] + 2 * gammaDotAlpha[i][j]) * dt3 / 6
                    + (alpha4[i][j] + alpha2DotGamma[i][j]
                        + 2 * alphaDotGammaDotAlpha[i][j]
                        + 3 * gammaDotAlpha2[i][j]
                        + gamma2[i][j]) * dt4 / 24 ;
        }
    }

    // Xi
    quad3x3Matrix xiMatrix;
    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < 3; ++j){
            xiMatrix[i][j] = I[i][j] * dt2 / 2 + alpha[i][j] * dt3 / 6 
                            + (3 * alpha2[i][j] + gamma[i][j]) * dt4 / 24
                            + (4 * alpha3[i][j] + 2 * alphaDotGamma[i][j] + 2 * gammaDotAlpha[i][j]) * dt5 / 120
                            + (5 * alpha4[i][j] + 3 * alpha2DotGamma[i][j]
                                + 4 * alphaDotGammaDotAlpha[i][j]
                                + 3 * gammaDotAlpha2[i][j]
                                + gamma2[i][j]) * dt6 / 720 ;
        }
    }
    Xi = matrixMultiply(xiMatrix, gVector);

    // Phi
    quad3x3Matrix phiMatrix;
    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < 3; ++j){
            phiMatrix[i][j] = I[i][j] * dt + alpha[i][j] * dt2 / 2 
                            + (alpha2[i][j] + gamma[i][j]) * dt3 / 6        
                            + (alpha3[i][j] + alphaDotGamma[i][j] + 2 * gammaDotAlpha[i][j]) * dt4 / 24
                            + (alpha4[i][j] + alpha2DotGamma[i][j] + 2 * alphaDotGammaDotAlpha[i][j]
                                + 3 * gammaDotAlpha2[i][j] + gamma2[i][j]) * dt5 / 120 ; 
        }
    }
    Phi = matrixMultiply(phiMatrix, gVector);

    // Store the double-precision versions
    A64 = convert3x3MatrixToDouble(A);
    B64 = convert3x3MatrixToDouble(B);
    C64 = convert3x3MatrixToDouble(C);
    D64 = convert3x3MatrixToDouble(D);
    Xi64  = convertThreeVectorToDouble(Xi);
    Phi64 = convertThreeVectorToDouble(Phi);

}

doubleThreeVector AISAbcdPropagator::CalculateNewPos(doubleThreeVector pos0, doubleThreeVector vel0, __float128 dt){
    doubleThreeVector pos1;
    pos1 = matrixAdd(matrixMultiply(A64, pos0), matrixMultiply(B64, vel0));
    pos1 = matrixAdd(pos1, Xi64);
    return pos1;
}

doubleThreeVector AISAbcdPropagator::CalculateNewVel(doubleThreeVector pos0, doubleThreeVector vel0, __float128 dt){
    doubleThreeVector vel1;
    vel1 = matrixAdd(matrixMultiply(C64, pos0), matrixMultiply(D64, vel0));
    vel1 = matrixAdd(vel1, Phi64);
    return vel1;
}