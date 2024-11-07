#include "AISKinematicPropagator.hh"

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

    __float128 t0 = wavePacket->GetTime();

    // convert t0 and t1 to double (static cast)
    double t0_double = (double)t0;
    double t1_double = (double)t1;

    // compute new phase space coordinates
    std::array<doubleThreeVector, 2> newPosVel = CalculateNewPhaseSpaceCoords(t0_double, t1_double, currentPos, currentVel);
    doubleThreeVector newPos = newPosVel[0];
    doubleThreeVector newVel = newPosVel[1];

    wavePacket->SetPosition(newPos);
    wavePacket->SetVelocity(newVel);

    std::array<double,2> phaseDoubleRes = CalculateNewPhaseDouble(currentPhaseDouble, currentPhaseDoubleErr, currentPos, currentVel, t0_double, t1_double);
    double newPhaseDouble    = phaseDoubleRes[0];
    double newPhaseDoubleErr = phaseDoubleRes[1];

    // if double phase error above 1e-3 radians, print a warning
    if(newPhaseDoubleErr > 1e-3)
    {
        std::cerr << "Warning: phase error is above 1e-3 radians! (" << newPhaseDoubleErr << " radians)" << std::endl;
    }

    wavePacket->SetPhaseDouble(newPhaseDouble);
    wavePacket->SetPhaseDoubleError(newPhaseDoubleErr);

    if(fAddEnergyPhase && wavePacket->GetState() == 1)
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

    __float128 t0 = wavePacket->GetTime();

    // convert t0 and t1 to double (static cast)
    double t0_double = (double)t0;
    double t1_double = (double)t1;

    // compute new phase space coordinates
    std::array<doubleThreeVector, 2> newPosVel = CalculateNewPhaseSpaceCoordsLinearized(t0_double, t1_double, currentPos, currentVel);
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
    double reltol = 0.0;
    double abstol = 1e-3;
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
                                                                                  const doubleThreeVector& pos0, const doubleThreeVector& vel0)
{   // check if t0 == t1 (can happen here)
    if(t0==t1)
    {
        return {pos0, vel0};
    }
    else
    {
    // define the ode system
    FuncLinearizedParams* linearizedParams = new FuncLinearizedParams{pos0, vel0, this};
    gsl_odeiv2_system sys = {funcLinearized, nullptr, 6, linearizedParams};
    // setup the driver
    double reltol = 0.0;
    double abstol = 1e-3;
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
                                                                                  const doubleThreeVector& pos0, const doubleThreeVector& vel0)
{
    // get the solutions to the linearized Hamilton's equations
    std::array<doubleThreeVector, 2> newCoords = CalculateNewPhaseSpaceCoordsLinearized(t0, t1, pos0, vel0);
    doubleThreeVector pos = newCoords[0];
    doubleThreeVector vel = newCoords[1];
    // get the gradients and hessians
    doubleThreeVector dUdx = this->get_dUdx(pos0,vel0);
    doubleThreeVector dUdp = get_dUdp(pos0,vel0);
    double3x3Matrix d2Udxdx = get_d2Udxdx(pos0,vel0);
    double3x3Matrix d2Udxdp = get_d2Udxdp(pos0,vel0);
    double3x3Matrix d2Udpdp = get_d2Udpdp(pos0,vel0);
    double3x3Matrix I = {{
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}
    }};
    // compute the time derivatives (given by the linearized Hamilton's equations)
    doubleThreeVector term1 = matrixAdd(vel0,dUdp);
    doubleThreeVector term2 = dotProduct(I,matrixAdd(vel, scalarMultiply(vel0,-1)));
    doubleThreeVector term3 = dotProduct(d2Udpdp,matrixAdd(scalarMultiply(vel,massSr87),scalarMultiply(vel0,-massSr87)));
    doubleThreeVector term4 = dotProduct(transpose(d2Udxdp),matrixAdd(pos,scalarMultiply(pos0,-1)));
    doubleThreeVector dotPos = matrixAdd(term1,matrixAdd(term2,matrixAdd(term3,term4)));

    term1 = scalarMultiply(dUdx,-1/massSr87);
    term2 = scalarMultiply(dotProduct(d2Udxdx,matrixAdd(pos,scalarMultiply(pos0,-1))),-1/massSr87);
    term3 = scalarMultiply(dotProduct(d2Udxdp,matrixAdd(vel,scalarMultiply(vel0,-1))),-1);
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
    doubleThreeVector pos0 = linearizedParams->pos0;
    doubleThreeVector vel0 = linearizedParams->vel0;

    // start by computing dU/dx and dU/dp and the hessians
    doubleThreeVector dUdx = propagator->get_dUdx(pos0, vel0);
    doubleThreeVector dUdp = propagator->get_dUdp(pos0, vel0);
    double3x3Matrix d2Udxdx = propagator->get_d2Udxdx(pos0, vel0);
    double3x3Matrix d2Udxdp = propagator->get_d2Udxdp(pos0, vel0);
    double3x3Matrix d2Udpdx = transpose(d2Udxdp);
    double3x3Matrix d2Udpdp = propagator->get_d2Udpdp(pos0, vel0);

    // get the shifted coordinates
    doubleThreeVector pos = {y[0], y[1], y[2]};
    doubleThreeVector vel = {y[3], y[4], y[5]};
    doubleThreeVector posMinusPos0 = matrixAdd(pos,scalarMultiply(pos0,-1));
    doubleThreeVector velMinusVel0 = matrixAdd(vel,scalarMultiply(vel0,-1));

    // get the terms in the linearized Hamilton's equations
    doubleThreeVector d2UdpdpDotVel = dotProduct(d2Udpdp,velMinusVel0);
    doubleThreeVector d2UdpdxDotPos = dotProduct(d2Udpdx,posMinusPos0);
    doubleThreeVector d2UdxdpDotVel = dotProduct(d2Udxdp,velMinusVel0);
    doubleThreeVector d2UdxdxDotPos = dotProduct(d2Udxdx,posMinusPos0);

    // define the f vector
    f[0] = vel[0] + dUdp[0] + d2UdpdpDotVel[0] * massSr87 + d2UdpdxDotPos[0];
    f[1] = vel[1] + dUdp[1] + d2UdpdpDotVel[1] * massSr87 + d2UdpdxDotPos[1];
    f[2] = vel[2] + dUdp[2] + d2UdpdpDotVel[2] * massSr87 + d2UdpdxDotPos[2];
    f[3] = -dUdx[0] / massSr87 - d2UdxdxDotPos[0] / massSr87 - d2UdxdpDotVel[0];
    f[4] = -dUdx[1] / massSr87 - d2UdxdxDotPos[1] / massSr87 - d2UdxdpDotVel[1];
    f[5] = -dUdx[2] / massSr87 - d2UdxdxDotPos[2] / massSr87 - d2UdxdpDotVel[2];

    return GSL_SUCCESS;
}

double AISKinematicPropagator::get_L(const double& t, void *params)
{   
    // get the Lagrangian parameters
    LagrangianParams *p = (LagrangianParams *) params;
    double t0 = p->t0;
    doubleThreeVector pos0 = p->pos0;
    doubleThreeVector vel0 = p->vel0;

    // compute the Lagrangian
    std::array<doubleThreeVector, 2> newPosVel = CalculateNewPhaseSpaceCoords(t0, t, pos0, vel0);
    doubleThreeVector newPos = newPosVel[0];
    doubleThreeVector newVel = newPosVel[1];
    double U = get_U(newPos, newVel);
    double T = 0.5 * (newVel[0] * newVel[0] + newVel[1] * newVel[1] + newVel[2] * newVel[2]); // * massSr87;
    return T - U;
}

double AISKinematicPropagator::get_L_wrapper(double t, void *params)
{
    LagrangianParams *lagrangianParams = static_cast<LagrangianParams*>(params);
    AISKinematicPropagator* propagator = lagrangianParams->propagator;
    return propagator->get_L(t,params);
}

std::array<double,2> AISKinematicPropagator::get_Scl(const double& t0, const double& t1, const doubleThreeVector& pos0, const doubleThreeVector& vel0)
{
    // compute the action by integrating the Lagrangian between t0 and t1
    double workspace_size = 1000; // max number of subintervals
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (workspace_size);

    double result, error;

    gsl_function integrand;
    integrand.function = &get_L_wrapper;
    integrand.params = new LagrangianParams{t0, pos0, vel0, this};

    double epsabs = 1e-9;
    double epsrel = 0.0;
    int key = 1;

    gsl_integration_qag (&integrand, t0, t1, epsabs, epsrel, workspace_size, 
                         key, w, &result, &error);

    gsl_integration_workspace_free (w);

    return {result, error};
}

std::array<double,2> AISKinematicPropagator::CalculateNewPhaseDouble(const double& phase0, const double& phaseErr, const doubleThreeVector& pos0, const doubleThreeVector& vel0,
                                                                          const double t0, const double t1)
{
    std::array<double,2> scl = get_Scl(t0, t1, pos0, vel0); // action and error
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