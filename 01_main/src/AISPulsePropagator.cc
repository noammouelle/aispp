#include "AISPulsePropagator.hh"

AISPulsePropagator::AISPulsePropagator(std::shared_ptr<AISLaserBeam> beam, 
                                       __float128 t0, __float128 t1,
                                       std::shared_ptr<AISKinematicPropagator> prop, double threshold)
    : laserBeam(beam), kinematicPropagator(prop), initTime(t0), finalTime(t1),
      initTimeDouble(static_cast<double>(t0)), finalTimeDouble(static_cast<double>(t1)), amplitudeThreshold(threshold) {
}

AISPulsePropagator::~AISPulsePropagator()
{
}

void AISPulsePropagator::PropagateEnsemble(std::unique_ptr<AISAtomEnsemble>& atomEnsemble)
{
    #pragma omp parallel for
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

        // print path and vertical velocity for debugging
        //std::cout<<wavepacket0->GetPath() << " " 
        //         << wavepacket0->GetVelocity()[2] << std::endl;

        // Set pos0 and vel0 and t0 (initial pos and vel of atom) for the new wavepacket
        wavepacket1->SetPos0(wavepacket0->GetPos0());
        wavepacket1->SetVel0(wavepacket0->GetVel0());
        wavepacket1->SetT0(wavepacket0->GetT0());

        // Set the position and velocity about which to linearize the operators
        wavepacket0->SetPosStar(wavepacket0->GetPosition());
        wavepacket0->SetVelStar(wavepacket0->GetVelocity());

        // Apply U1(t_0,t_0) transformation (do nothing)
        // Apply U2(t_0,t_0) transformation
        ApplyU2(wavepacket0, this->initTime);

        // Apply U3(t_0,t)
        if(this->useStaticApprox) // TODO: Remove. This is now obsolete.
        {
            ApplyU3StaticApprox(wavepacket0, wavepacket1, this->initTime, this->finalTime);
        }
        else if(this->ultraFast)
        {
            ApplyU3UltraFast(wavepacket0, wavepacket1, this->initTime, this->finalTime);
        }
        else
        {
            ApplyU3(wavepacket0, wavepacket1, this->initTime, this->finalTime);
        }

        // add the wavepackets to a temporary vector
        std::unique_ptr<wavePacketVector> newWavePacketsTemp(new wavePacketVector);
        newWavePacketsTemp->push_back(std::move(wavepacket0));
        newWavePacketsTemp->push_back(std::move(wavepacket1));

        // Apply the MC branching scheme if needed
        if(this->useMcBranching)
        {
            ApplyMCBranching(newWavePacketsTemp);
        }
        // Apply the path selection scheme if needed
        if(this->usePathSelection)
        {
            ApplyPathSelection(newWavePacketsTemp, pathsToSimulate);
        }
        // Apply the detectable volume selection scheme if needed
        if(this->useDetVolSelection)
        {
            ApplyDetVolSelection(newWavePacketsTemp);
        }
        
        // add the wavepackets to the vector if the amplitude is above the threshold
        ApplyCutoff(newWavePacketsTemp);

        // add the wavepackets to the new wavepacket vector
        for(int i = 0; i < newWavePacketsTemp->size(); ++i)
        {
            newWavePackets->push_back(std::move(newWavePacketsTemp->at(i)));
        }      
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

        // debugging output
        //std::cout<<"Path: " << currentWavepacket->GetPath() << " State: " 
        //             << currentWavepacket->GetState() << " Amplitude: " 
        //             << currentWavepacket->GetAmplitude() << "new Phase Double: "
        //                << currentWavepacket->GetPhaseDouble() << std::endl;
    }

    // debug
    //std::cout<<"Number of wavepackets left: "<<newWavePackets->size()<<std::endl;

    // Step 3: update the wavepacket vector
    atom->DeleteWavePackets();
    atom->AddWavePackets(newWavePackets);
}

void AISPulsePropagator::ApplyU1(std::unique_ptr<AISWavePacket>& wavepacket, __float128 t0, __float128 t1)
{
    // save the add energy phase flag and set to true
    bool addEnergyPhase = kinematicPropagator->fAddEnergyPhase;
    kinematicPropagator->SetAddEnergyPhase(true);

    // U1 is just kinematic propagation from t1 to t0 (note the reverse time order)
    wavepacket->SetTime(t1); // TODO: there must be a better way of doing this
    this->kinematicPropagator->PropagateWavePacket(wavepacket, t0);

    // restore the add energy phase flag
    kinematicPropagator->SetAddEnergyPhase(addEnergyPhase);
}

void AISPulsePropagator::ApplyU2(std::unique_ptr<AISWavePacket>& wavepacket, __float128 t0)
{
    if(wavepacket->GetState() == 1)
    {
        double currentPhaseDouble = wavepacket->GetPhaseDouble();
        __float128 currentPhaseQuad = wavepacket->GetPhaseQuad();
        
        doubleThreeVector pos = wavepacket->GetPosition();
        doubleThreeVector posStar = wavepacket->GetPosStar();
        doubleThreeVector kT0 = laserBeam->GetK(t0);
        __float128 omegaT0 = laserBeam->GetOmega(t0);
        double phi = laserBeam->GetPhi(posStar);
        doubleThreeVector gradPhi = laserBeam->GetDelPhi(posStar);

        // dot(kT0,pos) ≈ kz*z ~ 450 Mrad; keeping it in float64 causes ~65 µrad rounding noise
        // across the x0 grid. Store it in phaseQuad (__float128) to avoid precision loss.
        // dot(gradPhi,pos) is the wavefront-curvature term; it is small (~few rad) so float64 is
        // fine — keep it in phaseDouble so the MZI port phase relationship is preserved.
        double dPhaseDouble = - phi + dotProduct(matrixAdd(posStar, scalarMultiply(pos, -1.0)), gradPhi);
        __float128 dPhaseQuad = omegaT0 * t0 - static_cast<__float128>(dotProduct(kT0, pos));

        wavepacket->SetPhaseDouble(currentPhaseDouble + dPhaseDouble);
        wavepacket->SetPhaseQuad(currentPhaseQuad + dPhaseQuad);

        // shift the central momentum
        doubleThreeVector vel  = wavepacket->GetVelocity();
        doubleThreeVector dvel = scalarMultiply(matrixAdd(kT0,gradPhi),-hbar/massSr87);
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
        doubleThreeVector posStar = wavepacket->GetPosStar();
        doubleThreeVector velStar = wavepacket->GetVelStar();
        doubleThreeVector kT1 = laserBeam->GetK(t1);
        __float128 omegaT1 = laserBeam->GetOmega(t1);
        double phi = laserBeam->GetPhi(posStar);
        doubleThreeVector gradPhi = laserBeam->GetDelPhi(posStar);
        doubleThreeVector kPrime = matrixAdd(kT1,gradPhi);

        // Get the A and B matrices and the Xi vector
        std::tuple<double3x3Matrix, double3x3Matrix, doubleThreeVector> ABXi = kinematicPropagator->get_ABXi(t0,t1,posStar,velStar);
        double3x3Matrix A = std::get<0>(ABXi);
        double3x3Matrix B = std::get<1>(ABXi);
        doubleThreeVector Xi = std::get<2>(ABXi);

        double t0Double = static_cast<double>(t0);
        double t1Double = static_cast<double>(t1);

        doubleThreeVector AdotPos = dotProduct(A,pos);
        doubleThreeVector AdotKPrime = dotProduct(A,kPrime);
        AdotKPrime = scalarMultiply(AdotKPrime,-hbar/(2*massSr87));
        doubleThreeVector BdotAdotKPrime = dotProduct(B,AdotKPrime);
        double propTerm = dotProduct(kPrime,
                                     matrixAdd(matrixAdd(AdotPos, BdotAdotKPrime), Xi));

        // dot(kPrime, AdotPos) ≈ kz*z ~ 450 Mrad; keeping it in float64 causes ~65 µrad rounding
        // noise across the x0 grid. Store it in phaseQuad (__float128) to avoid precision loss.
        double dPhaseDouble = phi - dotProduct(posStar,gradPhi) + propTerm - dotProduct(kPrime, AdotPos);
        __float128 dPhaseQuad = - (omegaT1 - omegaSr87) * t1 - omegaSr87*t0
                                + static_cast<__float128>(dotProduct(kPrime, AdotPos));

        wavepacket->SetPhaseDouble(currentPhaseDouble + dPhaseDouble);
        wavepacket->SetPhaseQuad(currentPhaseQuad + dPhaseQuad);

        // shift the central phase space coordinates
        doubleThreeVector dPos = scalarMultiply(dotProduct(transpose(B), kPrime), -hbar/massSr87);
        doubleThreeVector dVel = scalarMultiply(dotProduct(transpose(A), kPrime), hbar/massSr87);
        
        wavepacket->SetVelocity(matrixAdd(vel,dVel));
        wavepacket->SetPosition(matrixAdd(pos,dPos));
    }
}

double AISPulsePropagator::getDelta(doubleThreeVector posPrime, doubleThreeVector velPrime, doubleThreeVector posStar, doubleThreeVector velStar, double t0, double t1, std::shared_ptr<AISLaserBeam> laserBeam)
{
    // get the effective wavevector
    doubleThreeVector kT1 = laserBeam->GetK(t1);
    doubleThreeVector kChirp = laserBeam->GetKChirp();
    doubleThreeVector gradPhi = laserBeam->GetDelPhi(posStar);
    doubleThreeVector kPrime = matrixAdd(kT1,gradPhi);

    // get the detuning
    __float128 omega   = laserBeam->GetOmega(0.0q);
    __float128 dOmega  = omega - omegaSr87;
    __float128 omegaChirp = laserBeam->GetFrequencyChirp();
    double dOmegaDouble = static_cast<double>(dOmega);

    // compute the recoil term
    double recoilTerm = dotProduct(kPrime,kPrime) * hbar / (2 * massSr87);

    // compute the time derivative of R hat '
    std::array<doubleThreeVector,2> dotCoordsPrime = kinematicPropagator->get_dotPhaseSpaceCoordsLinearized(t0,t1,posPrime,velPrime,posStar,velStar);
    doubleThreeVector posDotPrime = dotCoordsPrime[0];

    // compute the wavevector chirp term
    double kChirpTerm = 0.5 * dotProduct(kChirp,posPrime);

    // compute the frequency chirp term
    double omegaChirpTerm = - omegaChirp * t1;

    // compute the detuning
    double delta = dotProduct(kPrime,posDotPrime) - dOmegaDouble + recoilTerm + kChirpTerm + omegaChirpTerm;

    return delta;
}

double AISPulsePropagator::getDeltaUltraFast(doubleThreeVector posPrime, doubleThreeVector velPrime, doubleThreeVector posStar, doubleThreeVector velStar, double t0, double t1, std::shared_ptr<AISLaserBeam> laserBeam)
{
    // get the effective wavevector
    doubleThreeVector kT1 = laserBeam->GetK(t1);
    doubleThreeVector kChirp = laserBeam->GetKChirp();
    doubleThreeVector gradPhi = laserBeam->GetDelPhi(posStar);
    doubleThreeVector kPrime = matrixAdd(kT1,gradPhi);

    // get the detuning
    __float128 omega   = laserBeam->GetOmega(0.0q);
    __float128 dOmega  = omega - omegaSr87;
    __float128 omegaChirp = laserBeam->GetFrequencyChirp();
    double dOmegaDouble = static_cast<double>(dOmega);

    // compute the recoil term
    double recoilTerm = dotProduct(kPrime,kPrime) * hbar / (2 * massSr87);

    // compute the time derivative of R hat '
    doubleThreeVector posDotPrime = velPrime;

    // compute the wavevector chirp term
    double kChirpTerm = 0.5 * dotProduct(kChirp,posPrime);

    // compute the frequency chirp term
    double omegaChirpTerm = - omegaChirp * t1;

    // compute the detuning
    double delta = dotProduct(kPrime,posDotPrime) - dOmegaDouble + recoilTerm + kChirpTerm + omegaChirpTerm;

    return delta;
}

int AISPulsePropagator::funcU3(double t, const double y[], double f[], void *params)
{
    // extract the params
    U3Params* u3Params = static_cast<U3Params*>(params);
    doubleThreeVector pos0,vel0,posStar,velStar;
    std::shared_ptr<AISLaserBeam> laserBeam = u3Params->laserBeam;
    pos0 = u3Params->pos0;
    vel0 = u3Params->vel0;
    posStar = u3Params->posStar;
    velStar = u3Params->velStar;
    AISPulsePropagator* pulsePropagator = u3Params->pulsePropagator;
    double t0 = u3Params->t0;

    // Compute the linearly propagated phase space coordinates
    doubleThreeVector posPrime, velPrime;

    std::array<doubleThreeVector,2> newCoordsLinearized = pulsePropagator->kinematicPropagator->CalculateNewPhaseSpaceCoordsLinearized(t0,t,pos0,vel0,posStar,velStar);
    posPrime = newCoordsLinearized[0];
    velPrime = newCoordsLinearized[1];

    // compute the detuning 
    double delta = pulsePropagator->getDelta(posPrime,velPrime, posStar, velStar, t0,t,laserBeam);

    // compute the Rabi frequency
    double RabiFreq = laserBeam->GetRabiFreq(posPrime,t0,t);

     // start by computing dU/dx and dU/dp
    double reComplexAmplitdueGround = y[0];
    double reComplexAmplitdueExcited = y[1];
    double imComplexAmplitdueGround = y[2];
    double imComplexAmplitdueExcited = y[3];

    // detuning factor
    if (pulsePropagator->ignoreDetuning)
    {
        delta = 0.0;
    }

    // define the f vector
    f[0] =   RabiFreq / 2 * imComplexAmplitdueExcited;
    f[1] =   delta * imComplexAmplitdueExcited + RabiFreq / 2 * imComplexAmplitdueGround;
    f[2] = - RabiFreq / 2 * reComplexAmplitdueExcited;
    f[3] = - delta * reComplexAmplitdueExcited - RabiFreq / 2 * reComplexAmplitdueGround;

    return GSL_SUCCESS;

}

int AISPulsePropagator::funcU3StaticApprox(double t, const double y[], double f[], void *params)
{
    // extract the params
    U3Params* u3Params = static_cast<U3Params*>(params);
    doubleThreeVector pos0,vel0,posStar,velStar;
    std::shared_ptr<AISLaserBeam> laserBeam = u3Params->laserBeam;
    pos0 = u3Params->pos0;
    vel0 = u3Params->vel0;
    posStar = u3Params->posStar;
    velStar = u3Params->velStar;
    AISPulsePropagator* pulsePropagator = u3Params->pulsePropagator;
    double t0 = u3Params->t0;

    // Use static approximation: wavepacket does not move during the pulse (use the initial position and velocity instead of heisenberg picture)
    doubleThreeVector posPrime, velPrime;

    posPrime = pos0;
    velPrime = vel0;

    // compute the detuning 
    double delta = pulsePropagator->getDelta(posPrime,velPrime, posStar, velStar, t0,t,laserBeam);

    // compute the Rabi frequency
    double RabiFreq = laserBeam->GetRabiFreq(posPrime,t0,t);

     // start by computing dU/dx and dU/dp
    double reComplexAmplitdueGround = y[0];
    double reComplexAmplitdueExcited = y[1];
    double imComplexAmplitdueGround = y[2];
    double imComplexAmplitdueExcited = y[3];

    // detuning factor
    if (pulsePropagator->ignoreDetuning)
    {
        delta = 0.0;
    }

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
    U3Params* params = new U3Params{t0Double, 
                                    wavepacket0->GetPosition(), wavepacket0->GetVelocity(), 
                                    wavepacket0->GetPosStar(), wavepacket0->GetVelStar(),
                                    laserBeam, this};

    // define the ode system
    gsl_odeiv2_system sys = {funcU3, nullptr, 4, params};
    // setup the driver
    double reltol = this->relTol;
    double abstol = this->absTol;
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
        doubleThreeVector posStar = wavepacket0->GetPosStar();
        doubleThreeVector velStar = wavepacket0->GetVelStar();
        std::string path = wavepacket0->GetPath();
        std::vector<std::string> detectablePaths = wavepacket0->GetDetectablePaths();
        bool willInterfere = wavepacket0->GetWillInterfere();

        wavepacket0->SetAmplitude(amplitudeGround);
        wavepacket0->SetPhaseDouble(currentPhaseDouble+phaseGround);
        wavepacket0->SetPath(path + "0"); // append a 0 to the path
        wavepacket0->SetPosStar(posStar); // ensures the expansions are about the correct point for the next unitary transformation
        wavepacket0->SetVelStar(velStar);
        wavepacket0->SetDetectablePaths(detectablePaths);
        wavepacket0->SetWillInterfere(willInterfere);
  
        wavepacket1->SetAmplitude(amplitudeExcited);
        wavepacket1->SetPhaseDouble(currentPhaseDouble+phaseExcited);
        wavepacket1->SetPhaseQuad(currentPhaseQuad);
        wavepacket1->SetState(1);
        wavepacket1->SetPosition(pos);
        wavepacket1->SetVelocity(vel);
        wavepacket1->SetPath(path + "1"); // append a 1 to the path
        wavepacket1->SetPosStar(posStar); // ensures the expansions are about the correct point for the next unitary transformation
        wavepacket1->SetVelStar(velStar);
        wavepacket1->SetDetectablePaths(detectablePaths);
        wavepacket1->SetWillInterfere(willInterfere);
    }
    else
    {
        double currentPhaseDouble = wavepacket0->GetPhaseDouble();
        __float128 currentPhaseQuad = wavepacket0->GetPhaseQuad();
        doubleThreeVector pos = wavepacket0->GetPosition();
        doubleThreeVector vel = wavepacket0->GetVelocity();
        doubleThreeVector posStar = wavepacket0->GetPosStar();
        doubleThreeVector velStar = wavepacket0->GetVelStar();
        std::string path = wavepacket0->GetPath();
        std::vector<std::string> detectablePaths = wavepacket0->GetDetectablePaths();
        bool willInterfere = wavepacket0->GetWillInterfere();

        wavepacket0->SetAmplitude(amplitudeExcited);
        wavepacket0->SetPhaseDouble(currentPhaseDouble+phaseExcited);
        wavepacket0->SetPath(path + "1"); // append a 1 to the path
        wavepacket0->SetPosStar(posStar); // ensures the expansions are about the correct point for the next unitary transformation
        wavepacket0->SetVelStar(velStar);
        wavepacket0->SetDetectablePaths(detectablePaths);
        wavepacket0->SetWillInterfere(willInterfere);

        wavepacket1->SetAmplitude(amplitudeGround);
        wavepacket1->SetPhaseDouble(currentPhaseDouble+phaseGround);
        wavepacket1->SetPhaseQuad(currentPhaseQuad);
        wavepacket1->SetState(0);
        wavepacket1->SetPosition(pos);
        wavepacket1->SetVelocity(vel);
        wavepacket1->SetPath(path + "0"); // append a 0 to the path
        wavepacket1->SetPosStar(posStar); // ensures the expansions are about the correct point for the next unitary transformation
        wavepacket1->SetVelStar(velStar);
        wavepacket1->SetDetectablePaths(detectablePaths);
        wavepacket1->SetWillInterfere(willInterfere);
    }

    // free the params
    delete params;
}

void AISPulsePropagator::ApplyU3UltraFast(std::unique_ptr<AISWavePacket>& wavepacket0, std::unique_ptr<AISWavePacket>& wavepacket1, __float128 t0, __float128 t1)
{
    // define the params
    double t0Double = static_cast<double>(t0);
    double t1Double = static_cast<double>(t1);

    // Debug
    //std::cout << "Using ApplyU3UltraFast" << std::endl;

    // set the amplitudes and phases of the wavepackets
    if(wavepacket0->GetState()==0)
    {
        // get the wavepacket params
        double currentPhaseDouble = wavepacket0->GetPhaseDouble();
        __float128 currentPhaseQuad = wavepacket0->GetPhaseQuad();
        doubleThreeVector pos = wavepacket0->GetPosition();
        doubleThreeVector vel = wavepacket0->GetVelocity();
        doubleThreeVector posStar = wavepacket0->GetPosStar();
        doubleThreeVector velStar = wavepacket0->GetVelStar();
        std::string path = wavepacket0->GetPath();
        std::vector<std::string> detectablePaths = wavepacket0->GetDetectablePaths();
        bool willInterfere = wavepacket0->GetWillInterfere();
        double amplitudeGround = wavepacket0->GetAmplitude();

        // compute the detuning
        double detuning = 1*getDeltaUltraFast(pos, vel, posStar, velStar, t0, t1, laserBeam);
        // compute the Rabi frequency
        double Omega = laserBeam->GetRabiFreq(pos, t0, t1);
        // compute the trigonometric arg
        double trigArg = 0.5 * (t1Double - t0Double) * sqrt(Omega * Omega + detuning * detuning);
        // compute the complex exponential factor
        std::complex<double> expFactor = exp(std::complex<double>(0,-0.5 * (t1Double-t0Double) * detuning));

        // compute the complex amplitude of the g->e transition
        std::complex<double> Seg = expFactor * std::complex<double>(0, -1 * Omega / sqrt(Omega * Omega + detuning * detuning) * sin(trigArg));
        // compute the complex amplitude of the g->g transition
        std::complex<double> Sgg = expFactor * std::complex<double>(cos(trigArg), detuning / sqrt(Omega * Omega + detuning * detuning) * sin(trigArg));

        wavepacket0->SetAmplitude(amplitudeGround * abs(Sgg));
        wavepacket0->SetPhaseDouble(currentPhaseDouble + std::arg(Sgg));
        wavepacket0->SetPath(path + "0"); // append a 0 to the path
        wavepacket0->SetPosStar(posStar); // ensures the expansions are about the correct point for the next unitary transformation
        wavepacket0->SetVelStar(velStar);
        wavepacket0->SetDetectablePaths(detectablePaths);
        wavepacket0->SetWillInterfere(willInterfere);
  
        wavepacket1->SetAmplitude(amplitudeGround * abs(Seg));
        wavepacket1->SetPhaseDouble(currentPhaseDouble + std::arg(Seg));
        wavepacket1->SetPhaseQuad(currentPhaseQuad);
        wavepacket1->SetState(1);
        wavepacket1->SetPosition(pos);
        wavepacket1->SetVelocity(vel);
        wavepacket1->SetPath(path + "1"); // append a 1 to the path
        wavepacket1->SetPosStar(posStar); // ensures the expansions are about the correct point for the next unitary transformation
        wavepacket1->SetVelStar(velStar);
        wavepacket1->SetDetectablePaths(detectablePaths);
        wavepacket1->SetWillInterfere(willInterfere);

        // debugging output
        //std::cout << "path: " << path << " detuning: " << detuning << std::endl;

    }
    else
    {
        double currentPhaseDouble = wavepacket0->GetPhaseDouble();
        __float128 currentPhaseQuad = wavepacket0->GetPhaseQuad();
        doubleThreeVector pos = wavepacket0->GetPosition();
        doubleThreeVector vel = wavepacket0->GetVelocity();
        doubleThreeVector posStar = wavepacket0->GetPosStar();
        doubleThreeVector velStar = wavepacket0->GetVelStar();
        std::string path = wavepacket0->GetPath();
        std::vector<std::string> detectablePaths = wavepacket0->GetDetectablePaths();
        bool willInterfere = wavepacket0->GetWillInterfere();
        double amplitudeExcited = wavepacket0->GetAmplitude();

        // compute the detuning
        double detuning = 1*getDeltaUltraFast(pos, vel, posStar, velStar, t0, t1, laserBeam);
        // compute the Rabi frequency
        double Omega = laserBeam->GetRabiFreq(pos, t0, t1);
        // compute the trigonometric arg
        double trigArg = 0.5 * (t1Double - t0Double) * sqrt(Omega * Omega + detuning * detuning);
        // compute the complex exponential factor
        std::complex<double> expFactor = std::exp(std::complex<double>(0,-0.5*(t1Double-t0Double) * detuning));
        // compute the complex amplitude of the e->g transition
        std::complex<double> Sge = expFactor * std::complex<double>(0, -1 * Omega / sqrt(Omega * Omega + detuning * detuning) * sin(trigArg));
        // compute the complex amplitude of the e->e transition
        std::complex<double> See = expFactor * std::complex<double>(cos(trigArg), - 1 * detuning / sqrt(Omega * Omega + detuning * detuning) * sin(trigArg));

        wavepacket0->SetAmplitude(amplitudeExcited * abs(See));
        wavepacket0->SetPhaseDouble(currentPhaseDouble + std::arg(See));
        wavepacket0->SetPath(path + "1"); // append a 1 to the path
        wavepacket0->SetPosStar(posStar); // ensures the expansions are about the correct point for the next unitary transformation
        wavepacket0->SetVelStar(velStar);
        wavepacket0->SetDetectablePaths(detectablePaths);
        wavepacket0->SetWillInterfere(willInterfere);

        wavepacket1->SetAmplitude(amplitudeExcited * abs(Sge));
        wavepacket1->SetPhaseDouble(currentPhaseDouble + std::arg(Sge));
        wavepacket1->SetPhaseQuad(currentPhaseQuad);
        wavepacket1->SetState(0);
        wavepacket1->SetPosition(pos);
        wavepacket1->SetVelocity(vel);
        wavepacket1->SetPath(path + "0"); // append a 0 to the path
        wavepacket1->SetPosStar(posStar); // ensures the expansions are about the correct point for the next unitary transformation
        wavepacket1->SetVelStar(velStar);
        wavepacket1->SetDetectablePaths(detectablePaths);
        wavepacket1->SetWillInterfere(willInterfere);

        // debugging output
        //std::cout << "path: " << path << " detuning: " << detuning << std::endl;
    }
}

void AISPulsePropagator::ApplyU3StaticApprox(std::unique_ptr<AISWavePacket>& wavepacket0, std::unique_ptr<AISWavePacket>& wavepacket1, __float128 t0, __float128 t1)
{
    // define the params
    double t0Double = static_cast<double>(t0);
    double t1Double = static_cast<double>(t1);
    U3Params* params = new U3Params{t0Double, 
                                    wavepacket0->GetPosition(), wavepacket0->GetVelocity(), 
                                    wavepacket0->GetPosStar(), wavepacket0->GetVelStar(),
                                    laserBeam, this};

    // define the ode system
    gsl_odeiv2_system sys = {funcU3StaticApprox, nullptr, 4, params};
    // setup the driver
    double reltol = this->relTol;
    double abstol = this->absTol;
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
        doubleThreeVector posStar = wavepacket0->GetPosStar();
        doubleThreeVector velStar = wavepacket0->GetVelStar();
        std::string path = wavepacket0->GetPath();
        std::vector<std::string> detectablePaths = wavepacket0->GetDetectablePaths();
        bool willInterfere = wavepacket0->GetWillInterfere();

        wavepacket0->SetAmplitude(amplitudeGround);
        wavepacket0->SetPhaseDouble(currentPhaseDouble+phaseGround);
        wavepacket0->SetPath(path + "0"); // append a 0 to the path
        wavepacket0->SetPosStar(posStar); // ensures the expansions are about the correct point for the next unitary transformation
        wavepacket0->SetVelStar(velStar);
        wavepacket0->SetDetectablePaths(detectablePaths);
        wavepacket0->SetWillInterfere(willInterfere);
  
        wavepacket1->SetAmplitude(amplitudeExcited);
        wavepacket1->SetPhaseDouble(currentPhaseDouble+phaseExcited);
        wavepacket1->SetPhaseQuad(currentPhaseQuad);
        wavepacket1->SetState(1);
        wavepacket1->SetPosition(pos);
        wavepacket1->SetVelocity(vel);
        wavepacket1->SetPath(path + "1"); // append a 1 to the path
        wavepacket1->SetPosStar(posStar); // ensures the expansions are about the correct point for the next unitary transformation
        wavepacket1->SetVelStar(velStar);
        wavepacket1->SetDetectablePaths(detectablePaths);
        wavepacket1->SetWillInterfere(willInterfere);
    }
    else
    {
        double currentPhaseDouble = wavepacket0->GetPhaseDouble();
        __float128 currentPhaseQuad = wavepacket0->GetPhaseQuad();
        doubleThreeVector pos = wavepacket0->GetPosition();
        doubleThreeVector vel = wavepacket0->GetVelocity();
        doubleThreeVector posStar = wavepacket0->GetPosStar();
        doubleThreeVector velStar = wavepacket0->GetVelStar();
        std::string path = wavepacket0->GetPath();
        std::vector<std::string> detectablePaths = wavepacket0->GetDetectablePaths();
        bool willInterfere = wavepacket0->GetWillInterfere();

        wavepacket0->SetAmplitude(amplitudeExcited);
        wavepacket0->SetPhaseDouble(currentPhaseDouble+phaseExcited);
        wavepacket0->SetPath(path + "1"); // append a 1 to the path
        wavepacket0->SetPosStar(posStar); // ensures the expansions are about the correct point for the next unitary transformation
        wavepacket0->SetVelStar(velStar);
        wavepacket0->SetDetectablePaths(detectablePaths);
        wavepacket0->SetWillInterfere(willInterfere);

        wavepacket1->SetAmplitude(amplitudeGround);
        wavepacket1->SetPhaseDouble(currentPhaseDouble+phaseGround);
        wavepacket1->SetPhaseQuad(currentPhaseQuad);
        wavepacket1->SetState(0);
        wavepacket1->SetPosition(pos);
        wavepacket1->SetVelocity(vel);
        wavepacket1->SetPath(path + "0"); // append a 0 to the path
        wavepacket1->SetPosStar(posStar); // ensures the expansions are about the correct point for the next unitary transformation
        wavepacket1->SetVelStar(velStar);
        wavepacket1->SetDetectablePaths(detectablePaths);
        wavepacket1->SetWillInterfere(willInterfere);
    }

    // free the params
    delete params;
}


std::vector<std::string> AISPulsePropagator::GetInterferingPaths()
{
    return interferingPaths;
}
void AISPulsePropagator::SetInterferingPaths(std::vector<std::string> paths)
{
    interferingPaths = paths;
}

std::vector<std::string> AISPulsePropagator::GetPathsToSimulate()
{
    return pathsToSimulate;
}
void AISPulsePropagator::SetPathsToSimulate(std::vector<std::string> paths)
{
    pathsToSimulate = paths;
}

bool AISPulsePropagator::GetUsePathSelection()
{
    return usePathSelection;
}
void AISPulsePropagator::SetUsePathSelection(bool usePathSelectionValue)
{
    usePathSelection = usePathSelectionValue;
}

void AISPulsePropagator::SetIgnoreDetuning(bool ignoreDetuning)
{
    this->ignoreDetuning = ignoreDetuning;
}

void AISPulsePropagator::ApplyMCBranching(std::unique_ptr<wavePacketVector>& newWavePackets)
{
    // check if the wavepackets are children of an interfering wavepacket
    if (newWavePackets->at(0)->GetWillInterfere())
    {
        // check if the new wavepackets are on interfering paths, update willInterfere flag
        for (int i = 0; i < this->interferingPaths.size(); i++)
        {
            std::string interferingPath = this->interferingPaths[i];
            std::string wavepacket0Path = newWavePackets->at(0)->GetPath();
            std::string wavepacket1Path = newWavePackets->at(1)->GetPath();

            if (interferingPath.substr(0, wavepacket0Path.size()) == wavepacket0Path)
            {
            newWavePackets->at(0)->SetWillInterfere(true);
            break;
            }
            else
            {
            newWavePackets->at(0)->SetWillInterfere(false);
            }
        }

        for (int i = 0; i < this->interferingPaths.size(); i++)
        {
            std::string interferingPath = this->interferingPaths[i];
            std::string wavepacket1Path = newWavePackets->at(1)->GetPath();

            if (interferingPath.substr(0, wavepacket1Path.size()) == wavepacket1Path)
            {
            newWavePackets->at(1)->SetWillInterfere(true);
            break;
            }
            else
            {
            newWavePackets->at(1)->SetWillInterfere(false);
            }
        }
    }
    else
    {
        // the wavepackets are children of a non-interfering wavepacket
        // apply the MC branching scheme
        std::unique_ptr<AISWavePacket>& wavepacket0 = newWavePackets->at(0);
        std::unique_ptr<AISWavePacket>& wavepacket1 = newWavePackets->at(1);
        // compute the (non normalized) probability of each path
        double prob0_ = wavepacket0->GetAmplitude() * wavepacket0->GetAmplitude();
        double prob1_ = wavepacket1->GetAmplitude() * wavepacket1->GetAmplitude();
        // normalize the probabilities
        double prob0 = prob0_ / (prob0_ + prob1_);
        double prob1 = prob1_ / (prob0_ + prob1_);
        // sample which wavepacket to keep
        double r = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
        if (r < prob0)
        {
            // delete wavepacket 1
            newWavePackets->erase(newWavePackets->begin() + 1);
        }
        else
        {
            // delete wavepacket 0
            newWavePackets->erase(newWavePackets->begin());
        }
    }

}

void AISPulsePropagator::ApplyCutoff(std::unique_ptr<wavePacketVector>& newWavePackets)
{
    // create a new wavepacket vector
    std::unique_ptr<wavePacketVector> newWavePacketsTemp(new wavePacketVector);

    // loop over wavepackets in reverse order
    for(int i = 0; i < newWavePackets->size(); ++i)
    {
        // if the amplitude is above the threshold, add the wavepacket to the new vector
        std::unique_ptr<AISWavePacket>& currentWavepacket = newWavePackets->at(i);
        if(abs(currentWavepacket->GetAmplitude()) >= this->amplitudeThreshold)
        {
            newWavePacketsTemp->push_back(std::move(currentWavepacket));
        }
    }

    // update the wavepacket vector
    newWavePackets = std::move(newWavePacketsTemp);
}

void AISPulsePropagator::ApplyDetVolSelection(std::unique_ptr<wavePacketVector>& newWavePackets)
{
    // create a new wavepacket vector
    std::unique_ptr<wavePacketVector> newWavePacketsTemp(new wavePacketVector);

    // loop over wavepackets
    for(int i = 0; i < newWavePackets->size(); ++i)
    {
        // if the wavepacket is not on a detectable path, delete the wavepacket
        std::unique_ptr<AISWavePacket>& currentWavepacket = newWavePackets->at(i);
        std::vector<std::string> detectablePaths = currentWavepacket->GetDetectablePaths();
        std::string currentPath = currentWavepacket->GetPath();
        // loop over detectable paths
        bool isDetectable = false;
        for(int j = 0; j < detectablePaths.size(); ++j)
        {
            if(detectablePaths[j].substr(0,currentPath.size()) == currentPath)
            {
                isDetectable = true;
                break;
            }
        }
        // if the wavepacket is  detectable, keep the wavepacket
        if(isDetectable)
        {
            newWavePacketsTemp->push_back(std::move(currentWavepacket));
        }
    }

    // update the wavepacket vector
    newWavePackets = std::move(newWavePacketsTemp);
}

void AISPulsePropagator::ApplyPathSelection(std::unique_ptr<wavePacketVector>& newWavePackets, std::vector<std::string> pathsToSimulate)
{
    // create a new wavepacket vector
    std::unique_ptr<wavePacketVector> newWavePacketsTemp(new wavePacketVector);

    // loop over wavepackets
    for(int i = 0; i < newWavePackets->size(); ++i)
    {
        // if the wavepacket is on a path to simulate, keep the wavepacket
        std::unique_ptr<AISWavePacket>& currentWavepacket = newWavePackets->at(i);
        std::string currentPath = currentWavepacket->GetPath();
        // loop over paths to simulate
        bool isPathToSimulate = false;
        for(int j = 0; j < pathsToSimulate.size(); ++j)
        {
            if(pathsToSimulate[j].substr(0,currentPath.size()) == currentPath)
            {
                isPathToSimulate = true;
                break;
            }
        }
        // if the wavepacket is on a path to simulate, keep the wavepacket
        if(isPathToSimulate)
        {
            newWavePacketsTemp->push_back(std::move(currentWavepacket));
        }
    }

    // update the wavepacket vector
    newWavePackets = std::move(newWavePacketsTemp);
}

bool AISPulsePropagator::GetUseMcBranching()
{
    return useMcBranching;
}

void AISPulsePropagator::SetUseMcBranching(bool useMcBranchingVal)
{
    useMcBranching = useMcBranchingVal;
}

bool AISPulsePropagator::GetUseDetVolSelection()
{
    return useDetVolSelection;
}

void AISPulsePropagator::SetUseDetVolSelection(bool useDetVolSelectionVal)
{
    useDetVolSelection = useDetVolSelectionVal;
}

bool AISPulsePropagator::GetUseStaticApprox()
{
    return useStaticApprox;
}
void AISPulsePropagator::SetUseStaticApprox(bool useStaticApproxVal)
{
    useStaticApprox = useStaticApproxVal;
}