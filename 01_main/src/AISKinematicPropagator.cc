#include "AISKinematicPropagator.hh"

#ifdef USE_CUDA
#include "AISKinematicPropagatorGPU.hh"
#endif

AISKinematicPropagator::AISKinematicPropagator(std::shared_ptr<potentialFunctionType> aU, std::shared_ptr<gradPotentialFunctionType> aDUdx,
                                               std::shared_ptr<gradPotentialFunctionType> aDUdp,
                                               std::shared_ptr<hessianPotentialFunctionType> aD2Udxdx,
                                               std::shared_ptr<hessianPotentialFunctionType> aD2Udxdp,
                                               std::shared_ptr<hessianPotentialFunctionType> aD2Udpdp) : U(aU), dUdx(aDUdx), dUdp(aDUdp),
                                                                                                  d2Udxdx(aD2Udxdx), d2Udxdp(aD2Udxdp), d2Udpdp(aD2Udpdp)
{}

AISKinematicPropagator::~AISKinematicPropagator()
{}

double AISKinematicPropagator::get_U(const doubleThreeVector& pos, const doubleThreeVector& vel)
{
    if(U == nullptr)
    {
        std::cerr << "Potential function not initialized!" << std::endl;
        throw std::runtime_error("Potential function not initialized");
    }
    return (*U)(pos, vel);
}

doubleThreeVector AISKinematicPropagator::get_dUdx(const doubleThreeVector& pos, const doubleThreeVector& vel)
{
    if(dUdx == nullptr)
    {
        std::cerr << "Gradient potential function not initialized!" << std::endl;
        throw std::runtime_error("Gradient potential function not initialized");
    }
    return (*dUdx)(pos, vel);
}

doubleThreeVector AISKinematicPropagator::get_dUdp(const doubleThreeVector& pos, const doubleThreeVector& vel)
{
    if(dUdp == nullptr)
    {
        std::cerr << "Gradient potential function not initialized!" << std::endl;
        throw std::runtime_error("Gradient potential function not initialized");
    }
    return (*dUdp)(pos, vel);
}

double3x3Matrix AISKinematicPropagator::get_d2Udxdx(const doubleThreeVector& pos, const doubleThreeVector& vel)
{
    if(d2Udxdx == nullptr)
    {
        std::cerr << "Hessian potential function not initialized!" << std::endl;
        throw std::runtime_error("Hessian potential function not initialized");
    }
    return (*d2Udxdx)(pos, vel);
}

double3x3Matrix AISKinematicPropagator::get_d2Udxdp(const doubleThreeVector& pos, const doubleThreeVector& vel)
{
    if(d2Udxdp == nullptr)
    {
        std::cerr << "Hessian potential function not initialized!" << std::endl;
        throw std::runtime_error("Hessian potential function not initialized");
    }
    return (*d2Udxdp)(pos, vel);
}

double3x3Matrix AISKinematicPropagator::get_d2Udpdp(const doubleThreeVector& pos, const doubleThreeVector& vel)
{
    if(d2Udpdp == nullptr)
    {
        std::cerr << "Hessian potential function not initialized!" << std::endl;
        throw std::runtime_error("Hessian potential function not initialized");
    }
    return (*d2Udpdp)(pos, vel);
}

void AISKinematicPropagator::SetAddEnergyPhase(bool addEnergyPhase)
{
    fAddEnergyPhase = addEnergyPhase;
}

void AISKinematicPropagator::PropagateEnsemble(std::unique_ptr<AISAtomEnsemble>& atomEnsemble, __float128 t1)
{
#ifdef USE_CUDA
    // GPU path: ultrafast mode + constant-acceleration potential (zero or linear gravity).
    // For quadratic_pot (position-dependent force) we fall back to CPU — the force
    // is not constant, so the analytic update does not apply.
    if (ultraFast && (potentialType == "zero_pot" || potentialType == "linear_pot")) {
        PropagateEnsembleGPU(atomEnsemble, t1);
        return;
    }
#endif
    #pragma omp parallel for
    for(int i_atom = 0; i_atom < atomEnsemble->GetNumberOfAtoms(); ++i_atom)
    {
        std::unique_ptr<AISAtom>& currentAtom = atomEnsemble->GetAtom(i_atom);
        PropagateAtom(currentAtom, t1);
    }
}

void AISKinematicPropagator::PropagateAtom(std::unique_ptr<AISAtom>& atom, __float128 t1)
{
    for(int i_wavePacket = 0; i_wavePacket < atom->GetNumberOfWavePackets(); ++i_wavePacket)
    {
        std::unique_ptr<AISWavePacket>& currentWavePacket = atom->GetWavePacket(i_wavePacket);
        PropagateWavePacket(currentWavePacket, t1);
    }
}

void AISKinematicPropagator::PropagateWavePacket(std::unique_ptr<AISWavePacket>& wavePacket, __float128 t1)
{
    doubleThreeVector currentPos  = wavePacket->GetPosition();
    doubleThreeVector currentVel  = wavePacket->GetVelocity();

    double currentPhaseDouble    = wavePacket->GetPhaseDouble();
    double currentPhaseDoubleErr = wavePacket->GetPhaseDoubleError();

    // TODO: Remove this
    // apply mod 2pi to the phase double
    //currentPhaseDouble = fmod(currentPhaseDouble, 2 * M_PI);
    //wavePacket->SetPhaseDouble(currentPhaseDouble);

    // apply mod 2pi to the phase quad (quadruple precision)
    //__float128 currentPhaseQuad = wavePacket->GetPhaseQuad();
    //currentPhaseQuad = fmod(currentPhaseQuad, 2 * M_PI);
    //wavePacket->SetPhaseQuad(currentPhaseQuad);

    __float128 t0 = wavePacket->GetTime();

    doubleThreeVector atomInitialPos = wavePacket->GetPos0();
    doubleThreeVector atomInitialVel = wavePacket->GetVel0();
    __float128 atomInitialTime = wavePacket->GetT0();

    // convert t0 and t1 to double (static cast)
    double t0_double = (double)t0;
    double t1_double = (double)t1;

    // compute new phase space coordinates
    std::array<doubleThreeVector, 2> newPosVel = CalculateNewPhaseSpaceCoords(t0_double, t1_double, currentPos, currentVel);
    doubleThreeVector newPos = newPosVel[0];
    doubleThreeVector newVel = newPosVel[1];

    wavePacket->SetPosition(newPos);
    wavePacket->SetVelocity(newVel);

    std::array<double,2> phaseDoubleRes;

    if (ultraFast==false){
        phaseDoubleRes = CalculateNewPhaseDouble(currentPhaseDouble, currentPhaseDoubleErr, currentPos, currentVel,
                                                                  atomInitialPos, atomInitialVel,
                                                                  t0_double, t1_double, atomInitialTime);
    }
    else{
        phaseDoubleRes = {currentPhaseDouble, 0}; // ultra-fast version, ignore action phase
    }

    double newPhaseDouble    = phaseDoubleRes[0];
    double newPhaseDoubleErr = phaseDoubleRes[1];

    // if double phase error above 0.5e-3 radians, print a warning
    if(newPhaseDoubleErr > 5e-4)
    {
        std::cerr << "Warning: phase error is above 0.5e-3 radians! (" << newPhaseDoubleErr << " radians)" << std::endl;
    }

    wavePacket->SetPhaseDouble(newPhaseDouble);
    wavePacket->SetPhaseDoubleError(newPhaseDoubleErr);

    if(wavePacket->GetState() == 1)
    {
        wavePacket->SetPhaseQuad(CalculateNewPhaseQuad(wavePacket->GetPhaseQuad(), t0, t1));
    };
    
    wavePacket->SetTime(t1);
}

void AISKinematicPropagator::PropagateEnsembleLinearized(std::unique_ptr<AISAtomEnsemble>& atomEnsemble, __float128 t1)
{
    #pragma omp parallel for
    for(int i_atom = 0; i_atom < atomEnsemble->GetNumberOfAtoms(); ++i_atom)
    {
        std::unique_ptr<AISAtom>& currentAtom = atomEnsemble->GetAtom(i_atom);
        PropagateAtomLinearized(currentAtom, t1);
    }
}

void AISKinematicPropagator::PropagateAtomLinearized(std::unique_ptr<AISAtom>& atom, __float128 t1)
{
    for(int i_wavePacket = 0; i_wavePacket < atom->GetNumberOfWavePackets(); ++i_wavePacket)
    {
        std::unique_ptr<AISWavePacket>& currentWavePacket = atom->GetWavePacket(i_wavePacket);
        PropagateWavePacketLinearized(currentWavePacket, t1);
    }
}

void AISKinematicPropagator::PropagateWavePacketLinearized(std::unique_ptr<AISWavePacket>& wavePacket, __float128 t1)
{
    // only propagate kinematics, not the phase
    doubleThreeVector currentPos  = wavePacket->GetPosition();
    doubleThreeVector currentVel  = wavePacket->GetVelocity();
    doubleThreeVector currentPosStar = wavePacket->GetPosStar();
    doubleThreeVector currentVelStar = wavePacket->GetVelStar();

    __float128 t0 = wavePacket->GetTime();

    // convert t0 and t1 to double (static cast)
    double t0_double = (double)t0;
    double t1_double = (double)t1;

    // compute new phase space coordinates
    std::array<doubleThreeVector, 2> newPosVel = CalculateNewPhaseSpaceCoordsLinearized(t0_double, t1_double, currentPos, currentVel, currentPosStar, currentVelStar);
    doubleThreeVector newPos = newPosVel[0];
    doubleThreeVector newVel = newPosVel[1];

    wavePacket->SetPosition(newPos);
    wavePacket->SetVelocity(newVel);

    wavePacket->SetTime(t1);
}

std::array<doubleThreeVector, 2> AISKinematicPropagator::CalculateNewPhaseSpaceCoords(const double& t0, const double& t1, 
                                                                                  const doubleThreeVector& pos0, const doubleThreeVector& vel0)
{
    // define the ode system
    gsl_odeiv2_system sys = {func, nullptr, 6, this};
    // setup the driver
    double reltol = odeRelTol;
    double abstol = odeAbsTol;
    double hstart = (t1 - t0) / 10.0;
    gsl_odeiv2_driver * d =
    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
                                  hstart, abstol, reltol);
    double y[6] = {pos0[0], pos0[1], pos0[2], vel0[0], vel0[1], vel0[2]};
    
    double t_start = t0;
    int status = gsl_odeiv2_driver_apply(d, &t_start, t1, y);
    if (status != GSL_SUCCESS)
    {
        printf("error, return value=%d\n", status);
    }
    gsl_odeiv2_driver_free(d);

    doubleThreeVector newPos = {y[0], y[1], y[2]};
    doubleThreeVector newVel = {y[3], y[4], y[5]};

    return {newPos, newVel};
}

std::array<doubleThreeVector, 2> AISKinematicPropagator::CalculateNewPhaseSpaceCoordsLinearized(const double& t0, const double& t1, 
                                                                                  const doubleThreeVector& pos0, const doubleThreeVector& vel0,
                                                                                  const doubleThreeVector& posStar, const doubleThreeVector& velStar)
{   // check if t0 == t1 (can happen here)
    if(t0==t1)
    {
        return {pos0, vel0};
    }
    else
    {
    // define the ode system
    FuncLinearizedParams* linearizedParams = new FuncLinearizedParams{posStar, velStar, this};
    gsl_odeiv2_system sys = {funcLinearized, nullptr, 6, linearizedParams};
    // setup the driver
    double reltol = odeRelTol;
    double abstol = odeAbsTol;
    double hstart = (t1 - t0) / 10.0;   
    gsl_odeiv2_driver * d =
    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
                                  hstart, abstol, reltol);
    double y[6] = {pos0[0], pos0[1], pos0[2], vel0[0], vel0[1], vel0[2]};
    
    double t_start = t0;
    int status = gsl_odeiv2_driver_apply(d, &t_start, t1, y);
    if (status != GSL_SUCCESS)
    {
        printf("error, return value=%d\n", status);
    }
    gsl_odeiv2_driver_free(d);

    doubleThreeVector newPos = {y[0], y[1], y[2]};
    doubleThreeVector newVel = {y[3], y[4], y[5]};

    return {newPos, newVel};
    }
}

std::array<doubleThreeVector,2> AISKinematicPropagator::get_dotPhaseSpaceCoordsLinearized(const double& t0, const double& t1, 
                                                                                          const doubleThreeVector& posPrime, const doubleThreeVector& velPrime,
                                                                                          const doubleThreeVector& posStar, const doubleThreeVector& velStar)
{
    // get the gradients and hessians
    doubleThreeVector dUdx = this->get_dUdx(posStar,velStar);
    doubleThreeVector dUdp = get_dUdp(posStar,velStar);
    double3x3Matrix d2Udxdx = get_d2Udxdx(posStar,velStar);
    double3x3Matrix d2Udxdp = get_d2Udxdp(posStar,velStar);
    double3x3Matrix d2Udpdp = get_d2Udpdp(posStar,velStar);
    double3x3Matrix I = {{
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}
    }};
    // compute the time derivatives (given by the linearized Hamilton's equations)
    doubleThreeVector term1 = matrixAdd(velStar,dUdp);
    doubleThreeVector term2 = dotProduct(I,matrixAdd(velPrime, scalarMultiply(velStar,-1)));
    doubleThreeVector term3 = dotProduct(d2Udpdp,matrixAdd(scalarMultiply(velPrime,massSr87),scalarMultiply(velStar,-massSr87)));
    doubleThreeVector term4 = dotProduct(transpose(d2Udxdp),matrixAdd(posPrime,scalarMultiply(posStar,-1)));
    doubleThreeVector dotPos = matrixAdd(term1,matrixAdd(term2,matrixAdd(term3,term4)));

    term1 = scalarMultiply(dUdx,-1/massSr87);
    term2 = scalarMultiply(dotProduct(d2Udxdx,matrixAdd(posPrime,scalarMultiply(posStar,-1))),-1/massSr87);
    term3 = scalarMultiply(dotProduct(d2Udxdp,matrixAdd(velPrime,scalarMultiply(velStar,-1))),-1);
    doubleThreeVector dotVel = matrixAdd(term1,matrixAdd(term2,term3));

    return {dotPos,dotVel};
}                                                                                  

int AISKinematicPropagator::func(double t, const double y[], double f[], void *params)
{
    // Function to be used by the ODE solver
    // defined such that dy_i/dt = f_i(t, y_1, y_2, ..., y_n)

    (void)(t); /* avoid unused parameter warning */
    //(void)(params); /* avoid unused parameter warning */
    AISKinematicPropagator* propagator = static_cast<AISKinematicPropagator*>(params);

    // start by computing dU/dx and dU/dp
    doubleThreeVector pos = {y[0], y[1], y[2]};
    doubleThreeVector vel = {y[3], y[4], y[5]};
    doubleThreeVector dUdx = propagator->get_dUdx(pos, vel);
    doubleThreeVector dUdp = propagator->get_dUdp(pos, vel);

    // define the f vector
    f[0] = vel[0] + dUdp[0];
    f[1] = vel[1] + dUdp[1];
    f[2] = vel[2] + dUdp[2];
    f[3] = -dUdx[0] / massSr87;
    f[4] = -dUdx[1] / massSr87;
    f[5] = -dUdx[2] / massSr87;

    return GSL_SUCCESS;
}

int AISKinematicPropagator::funcLinearized(double t, const double y[], double f[], void *params)
{
    // Function to be used by the ODE solver
    // defined such that dy_i/dt = f_i(t, y_1, y_2, ..., y_n)
    (void)(t); /* avoid unused parameter warning */
    FuncLinearizedParams* linearizedParams = static_cast<FuncLinearizedParams*>(params);
    AISKinematicPropagator* propagator = linearizedParams->propagator;
    doubleThreeVector posStar = linearizedParams->posStar;
    doubleThreeVector velStar = linearizedParams->velStar;

    // start by computing dU/dx and dU/dp and the hessians
    doubleThreeVector dUdx = propagator->get_dUdx(posStar,velStar);
    doubleThreeVector dUdp = propagator->get_dUdp(posStar,velStar);
    double3x3Matrix d2Udxdx = propagator->get_d2Udxdx(posStar,velStar);
    double3x3Matrix d2Udxdp = propagator->get_d2Udxdp(posStar,velStar);
    double3x3Matrix d2Udpdx = transpose(d2Udxdp);
    double3x3Matrix d2Udpdp = propagator->get_d2Udpdp(posStar,velStar);

    // get the shifted coordinates
    doubleThreeVector pos = {y[0], y[1], y[2]};
    doubleThreeVector vel = {y[3], y[4], y[5]};
    doubleThreeVector posMinusPosStar = matrixAdd(pos,scalarMultiply(posStar,-1));
    doubleThreeVector velMinusVelStar = matrixAdd(vel,scalarMultiply(velStar,-1));

    // get the terms in the linearized Hamilton's equations
    doubleThreeVector d2UdpdpDotVel = dotProduct(d2Udpdp,velMinusVelStar);
    doubleThreeVector d2UdpdxDotPos = dotProduct(d2Udpdx,posMinusPosStar);
    doubleThreeVector d2UdxdpDotVel = dotProduct(d2Udxdp,velMinusVelStar);
    doubleThreeVector d2UdxdxDotPos = dotProduct(d2Udxdx,posMinusPosStar);

    // define the f vector
    f[0] = vel[0] + dUdp[0] + d2UdpdpDotVel[0] * massSr87 + d2UdpdxDotPos[0];
    f[1] = vel[1] + dUdp[1] + d2UdpdpDotVel[1] * massSr87 + d2UdpdxDotPos[1];
    f[2] = vel[2] + dUdp[2] + d2UdpdpDotVel[2] * massSr87 + d2UdpdxDotPos[2];
    f[3] = -dUdx[0] / massSr87 - d2UdxdxDotPos[0] / massSr87 - d2UdxdpDotVel[0];
    f[4] = -dUdx[1] / massSr87 - d2UdxdxDotPos[1] / massSr87 - d2UdxdpDotVel[1];
    f[5] = -dUdx[2] / massSr87 - d2UdxdxDotPos[2] / massSr87 - d2UdxdpDotVel[2];

    return GSL_SUCCESS;
}

double AISKinematicPropagator::get_dL(const double& t, void *params)
{   
    // get the Lagrangian parameters
    LagrangianParams *p = (LagrangianParams *) params;
    double t0 = p->t0;
    doubleThreeVector pos0 = p->pos0;
    doubleThreeVector vel0 = p->vel0;

    // get the initial kinematics of the atom
    doubleThreeVector atomInitialPos = p->atomInitialPos;
    doubleThreeVector atomInitialVel = p->atomInitialVel;
    __float128 atomInitialTime = p->atomInitialTime;
    double atomInitialTimeDouble = (double)atomInitialTime;

    // compute the exact trajectory of the wavepacket and the unperturbed
    // atom trajectory in a quadratic Hamiltonian
    std::array<doubleThreeVector, 2> newPosVel = CalculateNewPhaseSpaceCoords(t0, t, pos0, vel0);
    std::array<doubleThreeVector, 2> newPosVelTilde = CalculateNewPhaseSpaceCoordsLinearized(atomInitialTimeDouble, t, atomInitialPos, atomInitialVel, atomInitialPos, atomInitialVel);
    doubleThreeVector newPos = newPosVel[0];
    doubleThreeVector newVel = newPosVel[1];
    doubleThreeVector newPosTilde = newPosVelTilde[0];
    doubleThreeVector newVelTilde = newPosVelTilde[1];

    // compute the perturbations around the unperturbed atom trajectory in a quadratic Hamiltonian
    doubleThreeVector dx = matrixAdd(newPos,scalarMultiply(newPosTilde,-1));
    doubleThreeVector dv = matrixAdd(newVel,scalarMultiply(newVelTilde,-1));

    // compute the terms in the perturbations about unperturbed trajectory's Lagrangian
    double dxDotDUdx = dotProduct(dx,get_dUdx(newPosTilde,newVelTilde)); // first order perturbations
    double dvDotDUdp = dotProduct(dv,get_dUdp(newPosTilde,newVelTilde));
    double dvDotVTilde = dotProduct(dv,newVelTilde);
    double dvDotDv   = dotProduct(dv,dv);
    double dvDotD2UdvdvDotDv = dotProduct(dv,dotProduct(get_d2Udpdp(newPosTilde,newVelTilde),dv));
    double dxDotD2UdxdvDotDv = dotProduct(dx,dotProduct(get_d2Udxdp(newPosTilde,newVelTilde),dv));
    double dxDotD2Udxdx = dotProduct(dx,dotProduct(get_d2Udxdx(newPosTilde,newVelTilde),dx));

    // compute and sum the perturbations up to second order(divided by m)
    double dL = - dxDotDUdx/massSr87 + (dvDotVTilde - dvDotDUdp)
                + 0.5 * (dvDotDv - dvDotD2UdvdvDotDv) 
                - dxDotD2UdxdvDotDv/massSr87
                - 0.5 * dxDotD2Udxdx/massSr87;

    return dL;
}

double AISKinematicPropagator::get_dL_wrapper(double t, void *params)
{
    LagrangianParams *lagrangianParams = static_cast<LagrangianParams*>(params);
    AISKinematicPropagator* propagator = lagrangianParams->propagator;
    return propagator->get_dL(t,params);
}

std::array<double,2> AISKinematicPropagator::get_dScl(const double& t0, const double& t1, const doubleThreeVector& pos0, const doubleThreeVector& vel0,
                                                     const doubleThreeVector& atomInitialPos, const doubleThreeVector& atomInitialVel,
                                                     const __float128& atomInitialTime)
{
    // compute the action by integrating the Lagrangian between t0 and t1
    double workspace_size = 1000; // max number of subintervals
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (workspace_size);

    double result, error;

    gsl_function integrand;
    integrand.function = &get_dL_wrapper;
    integrand.params = new LagrangianParams{t0, pos0, vel0, atomInitialPos, atomInitialVel, atomInitialTime, this};

    double epsabs = qagAbsTol;
    double epsrel = qagRelTol;
    int key = 6;

    gsl_integration_qag (&integrand, t0, t1, epsabs, epsrel, workspace_size, 
                         key, w, &result, &error);

    gsl_integration_workspace_free (w);

    return {result, error};
}

std::tuple<double3x3Matrix, double3x3Matrix, doubleThreeVector> AISKinematicPropagator::get_ABXi(const double& t0, const double& t1, 
                                                                             const doubleThreeVector& posStar, const doubleThreeVector& velStar)
{
    // Compute the A, B and Xi matrices

    // get the hessians
    double3x3Matrix d2Udxdx = get_d2Udxdx(posStar,velStar);
    double3x3Matrix d2Udxdp = get_d2Udxdp(posStar,velStar);
    double3x3Matrix d2Udpdp = get_d2Udpdp(posStar,velStar);
    double3x3Matrix d2Udpdx = transpose(d2Udxdp);

    // construct the kinematic matrix
    double6x6Matrix M = {{
        {d2Udpdx[0][0], d2Udpdx[0][1], d2Udpdx[0][2], 1 + d2Udpdp[0][0], d2Udpdp[0][1], d2Udpdp[0][2]},
        {d2Udpdx[1][0], d2Udpdx[1][1], d2Udpdx[1][2], d2Udpdp[1][0], 1 + d2Udpdp[1][1], d2Udpdp[1][2]},
        {d2Udpdx[2][0], d2Udpdx[2][1], d2Udpdx[2][2], d2Udpdp[2][0], d2Udpdp[2][1], 1 + d2Udpdp[2][2]},
        {-d2Udxdx[0][0] / massSr87, -d2Udxdx[0][1] / massSr87, -d2Udxdx[0][2] / massSr87, -d2Udxdp[0][0], -d2Udxdp[0][1], -d2Udxdp[0][2]},
        {-d2Udxdx[1][0] / massSr87, -d2Udxdx[1][1] / massSr87, -d2Udxdx[1][2] / massSr87, -d2Udxdp[1][0], -d2Udxdp[1][1], -d2Udxdp[1][2]},
        {-d2Udxdx[2][0] / massSr87, -d2Udxdx[2][1] / massSr87, -d2Udxdx[2][2] / massSr87, -d2Udxdp[2][0], -d2Udxdp[2][1], -d2Udxdp[2][2]}
    }};

    // compute M * (t1 - t0)
    M = scalarMultiply(M, t1 - t0);

    // Compute the matrix exponential of M * (t1 - t0)
    double6x6Matrix expM = expMatrix(M);

    // Extract A, B, and Xi from the matrix exponential
    double3x3Matrix A = {{
        {expM[0][0], expM[0][1], expM[0][2]},
        {expM[1][0], expM[1][1], expM[1][2]},
        {expM[2][0], expM[2][1], expM[2][2]}
    }};

    double3x3Matrix B = {{
        {expM[0][3], expM[0][4], expM[0][5]},
        {expM[1][3], expM[1][4], expM[1][5]},
        {expM[2][3], expM[2][4], expM[2][5]}
    }};

    // Now get Xi by solving the linear system for r(t0)=v(t0)=0
    std::array<doubleThreeVector,2> XiPhi = CalculateNewPhaseSpaceCoordsLinearized(t0, t0, {0,0,0}, {0,0,0}, posStar, velStar);
    doubleThreeVector Xi = XiPhi[0];

    return {A, B, Xi};
}

std::array<double,2> AISKinematicPropagator::CalculateNewPhaseDouble(const double& phase0, const double& phaseErr, const doubleThreeVector& pos0, const doubleThreeVector& vel0,
                                                                     const doubleThreeVector& atomInitialPos, const doubleThreeVector& atomInitialVel,
                                                                     const double t0, const double t1, const __float128& atomInitialTime)
{
    std::array<double,2> scl = get_dScl(t0, t1, pos0, vel0, atomInitialPos, atomInitialVel, atomInitialTime); // action and error
    double action = scl[0];
    double error = scl[1];

    // divide the action by hbar to get the phase and return result
    double dphase = action * massSr87overHbar;/// hbar;
    double dphaseError = error * massSr87overHbar; // / hbar;

    return {phase0 + dphase, sqrt(dphaseError * dphaseError + phaseErr * phaseErr)};
}

__float128 AISKinematicPropagator::CalculateNewPhaseQuad(const __float128& phase0, const __float128& t0, const __float128 t1)
{
    return phase0 - omegaSr87 * (t1 - t0);
}