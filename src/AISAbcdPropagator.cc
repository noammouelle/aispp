#include "AISConstants.hh"
#include "AISAbcdPropagator.hh"


AISAbcdPropagator::AISAbcdPropagator(__float128 dt){
    deltaTime   = dt;
    //deltaTime64 = convertScalarToDouble(deltaTime);
    deltaTime64 = convertScalarToDouble(dt);

    // Identity
     I = {{ {{1.0, 0.0, 0.0}},
            {{0.0, 1.0, 0.0}},
            {{0.0, 0.0, 1.0}} }};
    
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
    double3x3Matrix gamma2 = dotProduct(gamma, gamma);

    double3x3Matrix alpha2 = dotProduct(alpha, alpha);
    double3x3Matrix alpha3 = dotProduct(alpha2, alpha);
    double3x3Matrix alpha4 = dotProduct(alpha3,alpha);

    double3x3Matrix alphaDotGamma = dotProduct(alpha, gamma);
    double3x3Matrix gammaDotAlpha = dotProduct(gamma, alpha);

    double3x3Matrix alpha2DotGamma = dotProduct(alpha2, gamma);
    double3x3Matrix gammaDotAlpha2 = dotProduct(gamma, alpha2);

    double3x3Matrix alphaDotGammaDotAlpha = dotProduct(alphaDotGamma, alpha);

    double dt  = deltaTime64;
    double dt2 = dt * dt;
    double dt3 = dt2 * dt;
    double dt4 = dt3 * dt;
    double dt5 = dt4 * dt;
    double dt6 = dt5 * dt;

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
    double3x3Matrix xiMatrix;
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
    Xi = dotProduct(xiMatrix, gVector);

    // Phi
    double3x3Matrix phiMatrix;
    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < 3; ++j){
            phiMatrix[i][j] = I[i][j] * dt + alpha[i][j] * dt2 / 2 
                            + (alpha2[i][j] + gamma[i][j]) * dt3 / 6        
                            + (alpha3[i][j] + alphaDotGamma[i][j] + 2 * gammaDotAlpha[i][j]) * dt4 / 24
                            + (alpha4[i][j] + alpha2DotGamma[i][j] + 2 * alphaDotGammaDotAlpha[i][j]
                                + 3 * gammaDotAlpha2[i][j] + gamma2[i][j]) * dt5 / 120 ; 
        }
    }
    Phi = dotProduct(phiMatrix, gVector);
}

doubleThreeVector AISAbcdPropagator::CalculateNewPos(const doubleThreeVector& pos0, const doubleThreeVector& vel0, __float128 dt){
    doubleThreeVector pos1;
    pos1 = matrixAdd(dotProduct(A, pos0), dotProduct(B, vel0));
    pos1 = matrixAdd(pos1, Xi);
    return pos1;
}

doubleThreeVector AISAbcdPropagator::CalculateNewVel(const doubleThreeVector& pos0, const doubleThreeVector& vel0, __float128 dt){
    doubleThreeVector vel1;
    vel1 = matrixAdd(dotProduct(C, pos0), dotProduct(D, vel0));
    vel1 = matrixAdd(vel1, Phi);
    return vel1;
}

double AISAbcdPropagator::CalculateNewPhaseDouble(const double& phase0, const doubleThreeVector& pos0, const doubleThreeVector& pos1, 
                                                  const doubleThreeVector& vel0, const doubleThreeVector& vel1){

    double phase1;
    // precompute matrix products
    double gVector2                    = dotProduct(gVector, gVector);
    doubleThreeVector gVectorDotGamma  = dotProduct(gVector, gamma);
    double gVectorDotGammaDotgVector   = dotProduct(gVectorDotGamma, gVector);
    doubleThreeVector gVectorDotGamma2 = dotProduct(gVectorDotGamma, gamma);
    double gVectorDotGamma2DotgVector  = dotProduct(gVectorDotGamma2, gVector);

    phase1 =  0.5 * dotProduct(pos1,vel1) - 0.5 * dotProduct(pos0,vel0);
    phase1 += 0.5 * dotProduct(Phi,pos1)  - 0.5 * dotProduct(Xi,vel1);
    phase1 += 0.5 * (gVector2 * pow(deltaTime64, 3)/6 
                        + gVectorDotGammaDotgVector  * pow(deltaTime64, 5)/120
                        + gVectorDotGamma2DotgVector * pow(deltaTime64, 7)/5040); //ignoring alpha
    phase1 *= massSr87 / hbar;

    return phase0 + phase1;
}