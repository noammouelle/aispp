#include "AISPulsePropagator.hh"

AISPulsePropagator::AISPulsePropagator(std::shared_ptr<AISLaserBeam> beam, 
                                       __float128 t0, __float128 t1,
                                       std::shared_ptr<AISKinematicPropagator> prop)
    : laserBeam(beam), kinematicPropagator(prop), initTime(t0), finalTime(t1),
      initTimeDouble(static_cast<double>(t0)), finalTimeDouble(static_cast<double>(t1)) {
}

AISPulsePropagator::~AISPulsePropagator()
{
}

void AISPulsePropagator::PropagateEnsemble(std::unique_ptr<AISAtomEnsemble>& atomEnsemble)
{
    //#pragma omp parallel for
    for(int i_atom = 0; i_atom < atomEnsemble->GetNumberOfAtoms(); ++i_atom)
    {
        std::unique_ptr<AISAtom>& currentAtom = atomEnsemble->GetAtom(i_atom);
        PropagateAtom(currentAtom);
    }
}

void AISPulsePropagator::PropagateAtom(std::unique_ptr<AISAtom>& atom)
{   
    // create a new wavepacket vector
    std::unique_ptr<wavePacketVector> newWavePackets(new wavePacketVector);

    // Step 1: Û3 * Û2 * Û1 (forward transformations)
    for(int wavePacketIndex = 0; wavePacketIndex < atom->GetNumberOfWavePackets(); ++wavePacketIndex)
    {
        // Get the pointer to the wavepackets
        std::unique_ptr<AISWavePacket>& wavepacket0 = atom->GetWavePacket(wavePacketIndex);
        std::unique_ptr<AISWavePacket> wavepacket1(new AISWavePacket());

        // Apply U1(t_0,t_0) transformation (do nothing)
        // Apply U2(t_0,t_0) transformation
        ApplyU2(wavepacket0, this->initTime);

        // Apply U3(t_0,t)
        ApplyU3(wavepacket0, wavepacket1, this->initTime, this->finalTime);
        
        // add the wavepackets to the vector
        newWavePackets->push_back(std::move(wavepacket0));
        newWavePackets->push_back(std::move(wavepacket1));        
    }

    // Step 2: perform the inverse transformations
    // Û1^dagger * Û_2^dagger
    //for(int wavePacketIndex = 0; wavePacketIndex < atom->GetNumberOfWavePackets(); ++wavePacketIndex)
    for(int wavePacketIndex = 0; wavePacketIndex < newWavePackets->size(); ++wavePacketIndex)
    {
        std::unique_ptr<AISWavePacket>& currentWavepacket = newWavePackets->at(wavePacketIndex);
        ApplyU2Dagger(currentWavepacket, this->initTime, this->finalTime);
        ApplyU1(currentWavepacket, this->finalTime, this->initTime); // note the reverse time order
        // make sure the wavepacket time is correct
        currentWavepacket->SetTime(this->finalTime);
    }

    // Step 3: update the wavepacket vector
    atom->DeleteWavePackets();
    atom->AddWavePackets(newWavePackets);
}

void AISPulsePropagator::ApplyU1(std::unique_ptr<AISWavePacket>& wavepacket, __float128 t0, __float128 t1)
{
    // U1 is just kinematic propagation from t1 to t0 (note the reverse time order)
    wavepacket->SetTime(t1); // TODO: there must be a better way of doing this
    this->kinematicPropagator->PropagateWavePacket(wavepacket, t0);
}

void AISPulsePropagator::ApplyU2(std::unique_ptr<AISWavePacket>& wavepacket, __float128 t0)
{
    if(wavepacket->GetState() == 1)
    {
        double currentPhaseDouble = wavepacket->GetPhaseDouble();
        __float128 currentPhaseQuad = wavepacket->GetPhaseQuad();
        
        doubleThreeVector pos = wavepacket->GetPosition();
        doubleThreeVector k = laserBeam->GetK();
        __float128 omega = laserBeam->GetOmega();
        double phi = laserBeam->GetPhi(pos);
        doubleThreeVector gradPhi = laserBeam->GetDelPhi(pos);

        double dPhaseDouble = - dotProduct(matrixAdd(k,gradPhi),pos) - phi;
        __float128 dPhaseQuad = omega * t0;

        wavepacket->SetPhaseDouble(currentPhaseDouble + dPhaseDouble);
        wavepacket->SetPhaseQuad(currentPhaseQuad + dPhaseQuad);

        // shift the central momentum
        doubleThreeVector vel  = wavepacket->GetVelocity();
        doubleThreeVector dvel = scalarMultiply(matrixAdd(k,gradPhi),-hbar/massSr87);
        wavepacket->SetVelocity(matrixAdd(vel,dvel));
    }
}

void AISPulsePropagator::ApplyU2Dagger(std::unique_ptr<AISWavePacket>& wavepacket, __float128 t0, __float128 t1)
{
    if(wavepacket->GetState() == 1)
    {
        double currentPhaseDouble = wavepacket->GetPhaseDouble();
        __float128 currentPhaseQuad = wavepacket->GetPhaseQuad();
        
        doubleThreeVector pos = wavepacket->GetPosition();
        doubleThreeVector vel = wavepacket->GetVelocity();
        doubleThreeVector k = laserBeam->GetK();
        __float128 omega = laserBeam->GetOmega();
        double phi = laserBeam->GetPhi(pos);
        doubleThreeVector gradPhi = laserBeam->GetDelPhi(pos);

        double t0Double = static_cast<double>(t0);
        double t1Double = static_cast<double>(t1);
        std::array<doubleThreeVector,2> newCoordsLinearized = kinematicPropagator->CalculateNewPhaseSpaceCoordsLinearized(t0Double,t1Double,pos,vel);
        doubleThreeVector posLin = newCoordsLinearized[0];
        doubleThreeVector velLin = newCoordsLinearized[1];

        double dPhaseDouble = phi + dotProduct(matrixAdd(k,gradPhi),posLin);
        __float128 dPhaseQuad = - (omega - omegaSr87) * t1 - omegaSr87*t0;

        wavepacket->SetPhaseDouble(currentPhaseDouble + dPhaseDouble);
        wavepacket->SetPhaseQuad(currentPhaseQuad + dPhaseQuad);

        // shift the central momentum
        wavepacket->SetVelocity(matrixAdd(vel,scalarMultiply(matrixAdd(k,gradPhi),hbar/massSr87)));
    }
}

double AISPulsePropagator::getDelta(doubleThreeVector pos0, doubleThreeVector vel0, double t0, double t1, std::shared_ptr<AISLaserBeam> laserBeam)
{
    // get the effective wavevector
    doubleThreeVector k = laserBeam->GetK();
    doubleThreeVector gradPhi = laserBeam->GetDelPhi(pos0);
    doubleThreeVector kPrime = matrixAdd(k,gradPhi);

    // get the detuning
    __float128 omega = laserBeam->GetOmega();
    __float128 dOmega = omega - omegaSr87;
    double dOmegaDouble = static_cast<double>(dOmega);

    // compute the recoil term
    double recoilTerm = dotProduct(kPrime,kPrime) * hbar / (2 * massSr87);

    // compute the time derivative of R hat '
    std::array<doubleThreeVector,2> dotCoordsPrime = kinematicPropagator->get_dotPhaseSpaceCoordsLinearized(t0,t1,pos0,vel0);
    doubleThreeVector posDotPrime = dotCoordsPrime[0];

    // compute the detuning
    double delta = dotProduct(kPrime,posDotPrime) - dOmegaDouble + recoilTerm;

    return delta;
}

int AISPulsePropagator::funcU3(double t, const double y[], double f[], void *params)
{
    // extract the params
    U3Params* u3Params = static_cast<U3Params*>(params);
    doubleThreeVector pos,vel;
    std::shared_ptr<AISLaserBeam> laserBeam = u3Params->laserBeam;
    pos = u3Params->pos;
    vel = u3Params->vel;
    AISPulsePropagator* pulsePropagator = u3Params->pulsePropagator;
    double t0 = u3Params->t0;

    // compute the detuning 
    double delta = pulsePropagator->getDelta(pos,vel,t0,t,laserBeam);

    // compute the Rabi frequency
    double RabiFreq = laserBeam->GetRabiFreq(pos,t0,t);

     // start by computing dU/dx and dU/dp
    double reComplexAmplitdueGround = y[0];
    double reComplexAmplitdueExcited = y[1];
    double imComplexAmplitdueGround = y[2];
    double imComplexAmplitdueExcited = y[3];

    // define the f vector
    f[0] =   RabiFreq / 2 * imComplexAmplitdueExcited;
    f[1] =   delta * imComplexAmplitdueExcited + RabiFreq / 2 * imComplexAmplitdueGround;
    f[2] = - RabiFreq / 2 * reComplexAmplitdueExcited;
    f[3] = - delta * reComplexAmplitdueExcited - RabiFreq / 2 * reComplexAmplitdueGround;

    return GSL_SUCCESS;

}

void AISPulsePropagator::ApplyU3(std::unique_ptr<AISWavePacket>& wavepacket0, std::unique_ptr<AISWavePacket>& wavepacket1, __float128 t0, __float128 t1)
{
    // define the params
    double t0Double = static_cast<double>(t0);
    double t1Double = static_cast<double>(t1);
    U3Params* params = new U3Params{t0Double, wavepacket0->GetPosition(), wavepacket0->GetVelocity(), laserBeam, this};
    // define the ode system
    gsl_odeiv2_system sys = {funcU3, nullptr, 4, params};
    // setup the driver
    double reltol = 0.0;
    double abstol = 1e-9;
    double hstart = (t1Double - t0Double) / 10.0;
    gsl_odeiv2_driver * d =
    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
                                  hstart, abstol, reltol);

    double y[4];
    if(wavepacket0->GetState()==0)
    {
        double amplitudeGround = wavepacket0->GetAmplitude();
        y[0] = amplitudeGround;
        y[1] = 0.0;
        y[2] = 0.0;
        y[3] = 0.0;
    }
    else
    {
        double amplitudeExcited = wavepacket0->GetAmplitude();
        y[0] = 0.0;
        y[1] = amplitudeExcited;
        y[2] = 0.0;
        y[3] = 0.0;
    }
    int status = gsl_odeiv2_driver_apply(d, &t0Double, t1Double, y);

    if (status != GSL_SUCCESS)
    {
        printf("error, return value=%d\n", status);
    }
    
    gsl_odeiv2_driver_free(d);

    // get the amplitudes 
    double reComplexAmplitdueGround = y[0];
    double reComplexAmplitdueExcited = y[1];
    double imComplexAmplitdueGround = y[2];
    double imComplexAmplitdueExcited = y[3];

    // separate into amplitude and phase
    double amplitudeGround = sqrt(reComplexAmplitdueGround * reComplexAmplitdueGround + imComplexAmplitdueGround * imComplexAmplitdueGround);
    double amplitudeExcited = sqrt(reComplexAmplitdueExcited * reComplexAmplitdueExcited + imComplexAmplitdueExcited * imComplexAmplitdueExcited);
    double phaseGround = atan2(imComplexAmplitdueGround,reComplexAmplitdueGround);
    double phaseExcited = atan2(imComplexAmplitdueExcited,reComplexAmplitdueExcited);

    // set the amplitudes and phases of the wavepackets
    if(wavepacket0->GetState()==0)
    {
        double currentPhaseDouble = wavepacket0->GetPhaseDouble();
        __float128 currentPhaseQuad = wavepacket0->GetPhaseQuad();
        doubleThreeVector pos = wavepacket0->GetPosition();
        doubleThreeVector vel = wavepacket0->GetVelocity();
        wavepacket0->SetAmplitude(amplitudeGround);
        wavepacket0->SetPhaseDouble(currentPhaseDouble+phaseGround);
  
        wavepacket1->SetAmplitude(amplitudeExcited);
        wavepacket1->SetPhaseDouble(currentPhaseDouble+phaseExcited);
        wavepacket1->SetPhaseQuad(currentPhaseQuad);
        wavepacket1->SetState(1);
        wavepacket1->SetPosition(pos);
        wavepacket1->SetVelocity(vel);
    }
    else
    {
        double currentPhaseDouble = wavepacket0->GetPhaseDouble();
        __float128 currentPhaseQuad = wavepacket0->GetPhaseQuad();
        doubleThreeVector pos = wavepacket0->GetPosition();
        doubleThreeVector vel = wavepacket0->GetVelocity();
        wavepacket0->SetAmplitude(amplitudeExcited);
        wavepacket0->SetPhaseDouble(currentPhaseDouble+phaseExcited);

        wavepacket1->SetAmplitude(amplitudeGround);
        wavepacket1->SetPhaseDouble(currentPhaseDouble+phaseGround);
        wavepacket1->SetPhaseQuad(currentPhaseQuad);
        wavepacket1->SetState(0);
        wavepacket1->SetPosition(pos);
        wavepacket1->SetVelocity(vel);
    }

    // free the params
    delete params;
}
