#include "AIStttPropagator.hh"
#include "AISLinearGravityPropagator.hh"
#include "AISUtilities.hh"
#include "AISConstants.hh"

AIStttPropagator::AIStttPropagator(__float128 deltaTime) : AISPulsePropagator(deltaTime)
{
    fDeltaTime = deltaTime;
    fDeltaTime64 = convertScalarToDouble(deltaTime);
    fFreePropagator = new AISLinearGravityPropagator(deltaTime / 2);

    fOmega0 = omegaSr87;

    // set the hamiltonian coefficients using the free 
    // propagator
    alpha = fFreePropagator->alpha;
    gamma = fFreePropagator->gamma;
    gVector = fFreePropagator->gVector;

    // identity matrix
     I = {{ {{1.0, 0.0, 0.0}},
            {{0.0, 1.0, 0.0}},
            {{0.0, 0.0, 1.0}} }};

    // set the ABCD matrices
    A = computeA(alpha, gamma, gVector, fDeltaTime64 / 2);
    B = computeB(alpha, gamma, gVector, fDeltaTime64 / 2);
    C = computeC(alpha, gamma, gVector, fDeltaTime64 / 2);
    D = computeD(alpha, gamma, gVector, fDeltaTime64 / 2);
    Xi = computeXi(alpha, gamma, gVector, fDeltaTime64 / 2);
    Phi = computePhi(alpha, gamma, gVector, fDeltaTime64 / 2);

    A_ = computeA(alpha, gamma, gVector, -fDeltaTime64 / 2);
    B_ = computeB(alpha, gamma, gVector, -fDeltaTime64 / 2);
    C_ = computeC(alpha, gamma, gVector, -fDeltaTime64 / 2);
    D_ = computeD(alpha, gamma, gVector, -fDeltaTime64 / 2);
    Xi_ = computeXi(alpha, gamma, gVector, -fDeltaTime64 / 2);
    Phi_ = computePhi(alpha, gamma, gVector, -fDeltaTime64 / 2);
}

AIStttPropagator::~AIStttPropagator()
{
    delete fFreePropagator;
}

void AIStttPropagator::setAmplitudeThreshold(double A)
{
    amplitudeThreshold = A;
}

double AIStttPropagator::computeDetuning(const doubleThreeVector& pos, const doubleThreeVector& vel,
                                         const __float128& omega)
{
    double detuning; 

    double domega = convertScalarToDouble(omega - fOmega0);
    
    doubleThreeVector v1 = dotProduct(C_, pos);
    doubleThreeVector v2 = dotProduct(D_, vel);
    doubleThreeVector v3 = dotProduct(alpha, dotProduct(A_, pos));
    doubleThreeVector v4 = dotProduct(alpha, dotProduct(B_, vel));
    doubleThreeVector v5 = dotProduct(alpha, Xi_);

    doubleThreeVector v  = matrixAdd(v1, v2);
    v = matrixAdd(v, v3);
    v = matrixAdd(v, v4);
    v = matrixAdd(v, v5);
    v = matrixAdd(v, Phi_);
    v = matrixAdd(v, scalarMultiply(fK, hbar/(2*massSr87)));

    detuning = domega - dotProduct(fK, v);

    return detuning;
}

complexDouble AIStttPropagator::Smatrix(const doubleThreeVector& pos, const doubleThreeVector& vel,
                                 const __float128& omega, const double& detuning, 
                                 const std::string& transition)
{   
    // compute the detuning
    double computedDetuning = detuning;
    // compute the effectrive Rabi frequency
    double effectiveRabiFreq = fIntensityProfile->GetEffectiveRabiFreq(pos);
    // compute the S matrix elements
    double X   = sqrt(computedDetuning * computedDetuning + 4*effectiveRabiFreq*effectiveRabiFreq);
    double arg = X/4 * fDeltaTime64;

    complexDouble Sij;

    if(transition=="g->g"){
        Sij = complexDouble(cos(arg), -computedDetuning*sin(arg)/X);
    }
    if(transition=="e->e"){
        Sij = complexDouble(cos(arg), computedDetuning*sin(arg)/X);
    }
    if(transition=="g->e" or transition=="e->g"){
        Sij = complexDouble(0.0, 2*effectiveRabiFreq*sin(arg)/X);
    }

    return Sij;
}

void AIStttPropagator::PropagateAtom(AISAtom* atom)
{
    // propagate freely half a time step
    fFreePropagator->PropagateAtom(atom);
    // propagate the atom using the TTT propagator
    TttPropagateAtom(atom);
    // propagator freely another time step
    fFreePropagator->PropagateAtom(atom);
}

void AIStttPropagator::TttPropagateAtom(AISAtom* atom)
{
    // create an empty wavepacket vector
    wavePacketVector* newWavePackets = new wavePacketVector;

    // some dummy variables
    doubleThreeVector dx;
    doubleThreeVector dv;

    doubleThreeVector newPos1;
    doubleThreeVector newVel1;
    doubleThreeVector newPos2;
    doubleThreeVector newVel2;

    double newAmplitude1;
    double newAmplitude2;
    double newPhaseDouble1;
    double newPhaseDouble2;
    __float128 newPhaseQuad1;
    __float128 newPhaseQuad2;

    // loop over all the wavepackets in the atom
    for(int wavePacketIndex = 0; wavePacketIndex < atom->GetNumberOfWavePackets(); ++wavePacketIndex)
    {
        // Get pointers to the two output wavepackets
        AISWavePacket* currentWavePacket = atom->GetWavePacket(wavePacketIndex);
        AISWavePacket* wavePacket1 = new AISWavePacket();
        AISWavePacket* wavePacket2 = new AISWavePacket();

        // copy info from original waveacket
        wavePacket1->SetPosition(currentWavePacket->GetPosition());
        wavePacket1->SetVelocity(currentWavePacket->GetVelocity());
        wavePacket1->SetTime(currentWavePacket->GetTime());
        wavePacket1->SetAmplitude(currentWavePacket->GetAmplitude());
        wavePacket1->SetPhaseDouble(currentWavePacket->GetPhaseDouble());
        wavePacket1->SetPhaseQuad(currentWavePacket->GetPhaseQuad());
        wavePacket1->SetState(currentWavePacket->GetState());

        wavePacket2->SetTime(currentWavePacket->GetTime());

        // dummy vector
        doubleThreeVector v;
        // dummy scalars
        double term1, term2, term2a, term2b, term3, term4;
        __float128 quadTerm;

        // for ground wavepacket
        if(wavePacket1->GetState()==0){
            // wp2 is in the excited state
            wavePacket2->SetState(1);

            /* ground -> ground transition */
            // rotate the amplitude and phase
            double detuning = computeDetuning(wavePacket1->GetPosition(), wavePacket1->GetVelocity(), fOmega);
            complexDouble Sgg = Smatrix(wavePacket1->GetPosition(), wavePacket1->GetVelocity(), fOmega, detuning, "g->g");

            newPos1 = wavePacket1->GetPosition();
            newVel1 = wavePacket1->GetVelocity();

            newAmplitude1 = wavePacket1->GetAmplitude() * std::abs(Sgg);
            newPhaseDouble1 = wavePacket1->GetPhaseDouble() + std::arg(Sgg) + detuning/2 * fDeltaTime64;
            
            newPhaseQuad1 = wavePacket1->GetPhaseQuad();

            /* ground -> excited transition */
            // propagate the kinematics
            dx = scalarMultiply(dotProduct(transpose(B),fK), - hbar/massSr87);
            dv = scalarMultiply(dotProduct(transpose(A),fK), hbar/massSr87);

            newPos2 = matrixAdd(wavePacket1->GetPosition(), dx);
            newVel2 = matrixAdd(wavePacket1->GetVelocity(), dv);

            // rotate the amplitude and phase
            complexDouble Sge = Smatrix(wavePacket1->GetPosition(), wavePacket1->GetVelocity(), fOmega, detuning, "g->e");

            newAmplitude2 = wavePacket1->GetAmplitude() * std::abs(Sge);
            newPhaseDouble2 = wavePacket1->GetPhaseDouble() + std::arg(Sge) + detuning/2 * fDeltaTime64;

            // add the contribution of the EM field to the phase
            v = dotProduct(transpose(B),fK);
            v = dotProduct(A,v);
            term1 = - hbar / (2*massSr87) * dotProduct(fK, v);
            v = dotProduct(A, wavePacket1->GetPosition());
            term2 = dotProduct(fK, v);
            term3 = dotProduct(fK, Xi);
            term4 = - fWaveFront->GetValue(wavePacket1->GetPosition());
            quadTerm = - fOmega * (wavePacket1->GetTime() + fDeltaTime / 2);
            // increment the phase
            newPhaseDouble2 += term1 + term2 + term3 + term4;
            newPhaseQuad2 = wavePacket1->GetPhaseQuad() + quadTerm;
        }
        if(wavePacket1->GetState()==1){
            // wp2 is in the ground state
            wavePacket2->SetState(0);

            /* excited -> ground transition */
            // propagate the kinematics
            dx = scalarMultiply(dotProduct(transpose(B_),fK), hbar/massSr87);
            dv = scalarMultiply(dotProduct(transpose(A_),fK), - hbar/massSr87);

            newPos2 = matrixAdd(wavePacket1->GetPosition(), dx);
            newVel2 = matrixAdd(wavePacket1->GetVelocity(), dv);

            // rotate the amplitude and phase
            doubleThreeVector effectivePos = matrixAdd(wavePacket1->GetPosition(), 
                                                       scalarMultiply(dotProduct(transpose(B_), fK), hbar/massSr87));
            doubleThreeVector effectiveVel = matrixAdd(wavePacket1->GetVelocity(), 
                                                       scalarMultiply(dotProduct(transpose(A_), fK), -hbar/massSr87));

            double detuning = computeDetuning(effectivePos, effectiveVel, fOmega);
            complexDouble Seg = Smatrix(effectivePos, effectiveVel, fOmega, detuning, "e->g");

            newAmplitude2 = wavePacket1->GetAmplitude() * std::abs(Seg);
            newPhaseDouble2 = wavePacket1->GetPhaseDouble() + std::arg(Seg) + detuning/2 * fDeltaTime64;

            // add the contribution of the EM field to the phase
            //doubleThreeVector v;
            v = dotProduct(transpose(B_),fK);
            v = dotProduct(A_,v);
            term1 = - hbar / (2*massSr87) * dotProduct(fK, v);
            v = dotProduct(A_, wavePacket1->GetPosition());
            term2 = - dotProduct(fK, v);
            term3 = - dotProduct(fK, Xi_);
            term4 = fWaveFront->GetValue(wavePacket1->GetPosition());
            quadTerm = fOmega * (wavePacket1->GetTime() - fDeltaTime / 2);
            // increment the phase
            newPhaseDouble2 += term1 + term2 + term3 + term4;
            newPhaseQuad2 = wavePacket1->GetPhaseQuad() + quadTerm;

            long double quadTermLd = static_cast<long double>(quadTerm);
            char quadTermLdStr[50];  // Adjust the buffer size as needed
            snprintf(quadTermLdStr, sizeof(quadTermLdStr), "%.34Le", quadTermLd);

            /* excited -> excited transition */
            // propagate the kinematics
            dx = matrixAdd(scalarMultiply(dotProduct(transpose(B_), fK), hbar/massSr87),
                           scalarMultiply(dotProduct(transpose(B), fK), -hbar/massSr87));
            dv = matrixAdd(scalarMultiply(dotProduct(transpose(A), fK), hbar/massSr87),
                           scalarMultiply(dotProduct(transpose(A_),fK), -hbar/massSr87));

            newPos1 = matrixAdd(wavePacket1->GetPosition(), dx);
            newVel1 = matrixAdd(wavePacket1->GetVelocity(), dv);

            // rotate the amplitude and phase
            complexDouble See = Smatrix(effectivePos, effectiveVel, fOmega, detuning, "e->e");

            newAmplitude1 = wavePacket1->GetAmplitude() * std::abs(See);
            newPhaseDouble1 = wavePacket1->GetPhaseDouble() + std::arg(See) + detuning/2 * fDeltaTime64;
            
            // add the contribution of the EM field to the phase
            //doubleThreeVector v;
            v = dotProduct(matrixAdd(dotProduct(A, transpose(B)), dotProduct(A, transpose(B))),fK);
            term1 = - hbar / (2*massSr87) * dotProduct(fK,v);
            v = dotProduct(matrixAdd(A, scalarMultiply(A_, -1)), wavePacket1->GetPosition());
            term2a = dotProduct(fK, v);
            v = dotProduct(transpose(B_),fK);
            v = dotProduct(A,v);
            term2b = hbar / massSr87 * dotProduct(fK, v);
            term2 = term2a + term2b;
            v = matrixAdd(Xi, scalarMultiply(Xi_, -1));
            term3 = dotProduct(fK, v);
            quadTerm = - fOmega * fDeltaTime;

            // increment the phase
            newPhaseDouble1 += term1 + term2 + term3;
            newPhaseQuad1 = wavePacket1->GetPhaseQuad() + quadTerm;
        }

        // update the wavepacket parameters
        wavePacket1->SetPosition(newPos1);
        wavePacket1->SetVelocity(newVel1);
        wavePacket1->SetAmplitude(newAmplitude1);
        wavePacket1->SetPhaseDouble(newPhaseDouble1);
        wavePacket1->SetPhaseQuad(newPhaseQuad1);
        wavePacket2->SetPosition(newPos2);
        wavePacket2->SetVelocity(newVel2);
        wavePacket2->SetAmplitude(newAmplitude2);
        wavePacket2->SetPhaseDouble(newPhaseDouble2);
        wavePacket2->SetPhaseQuad(newPhaseQuad2);

        // add the new wavepackets to the vector if the amplitude is above the threshold
        if(newAmplitude1 > amplitudeThreshold){
            newWavePackets->push_back(wavePacket1);
        }
        else{
            delete wavePacket1;
        }
        if(newAmplitude2 > amplitudeThreshold){
            newWavePackets->push_back(wavePacket2);
        }
        else{
            delete wavePacket2;
        }
    // delete the old wavepackets
    atom->DeleteWavePackets();
    // add the new wavepackets to the atom
    atom->AddWavePackets(newWavePackets);
}