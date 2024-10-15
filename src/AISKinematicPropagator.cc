#include "AISKinematicPropagator.hh"

AISKinematicPropagator::AISKinematicPropagator()
{
}

AISKinematicPropagator::~AISKinematicPropagator()
{}

void AISKinematicPropagator::SetAddEnergyPhase(bool addEnergyPhase)
{
    fAddEnergyPhase = addEnergyPhase;
}

void AISKinematicPropagator::PropagateEnsemble(AISAtomEnsemble* atomEnsemble, __float128 t1)
{
    #pragma omp parallel for
    for(int i_atom = 0; i_atom < atomEnsemble->GetNumberOfAtoms(); ++i_atom)
    {
        AISAtom* currentAtom = atomEnsemble->GetAtom(i_atom);
        PropagateAtom(currentAtom, t1);
    }
}

void AISKinematicPropagator::PropagateAtom(AISAtom* atom, __float128 t1)
{
    for(int i_wavePacket = 0; i_wavePacket < atom->GetNumberOfWavePackets(); ++i_wavePacket)
    {
        AISWavePacket* currentWavePacket = atom->GetWavePacket(i_wavePacket);
        PropagateWavePacket(currentWavePacket, t1);
    }
}

void AISKinematicPropagator::PropagateWavePacket(AISWavePacket* wavePacket, __float128 t1)
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
    std::array<doubleThreeVector, 2> newPosVel = CalculateNewPhaseSpaceCoords(t0_double, t1_double, currentPos, currentPos, currentVel, currentVel);
    doubleThreeVector newPos = newPosVel[0];
    doubleThreeVector newVel = newPosVel[1];

    wavePacket->SetPosition(newPos);
    wavePacket->SetVelocity(newVel);

    std::array<double,2> phaseDoubleRes = CalculateNewPhaseDouble(currentPhaseDouble, currentPhaseDoubleErr, currentPos, currentVel, t0_double, t1_double);
    double newPhaseDouble    = phaseDoubleRes[0];
    double newPhaseDoubleErr = phaseDoubleRes[1];

    wavePacket->SetPhaseDouble(newPhaseDouble);
    wavePacket->SetPhaseDoubleError(newPhaseDoubleErr);

    if(fAddEnergyPhase && wavePacket->GetState() == 1)
    {
        wavePacket->SetPhaseQuad(CalculateNewPhaseQuad(wavePacket->GetPhaseQuad(), t0, t1));
    };
    
    wavePacket->SetTime(t1);
}

std::array<doubleThreeVector, 2> AISKinematicPropagator::CalculateNewPhaseSpaceCoords(const double& t0, const double& t1, 
                                                                                  const doubleThreeVector& pos0, const doubleThreeVector& pos1, 
                                                                                  const doubleThreeVector& vel0, const doubleThreeVector& vel1)
{
    // define the ode system
    gsl_odeiv2_system sys = {func, nullptr, 6, this};
    // setup the driver
    double reltol = 0.0;
    double abstol = 1e-9;
    double hstart = (t1 - t0) / 1000.0;
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

double AISKinematicPropagator::get_L(const double& t, void *params)
{   
    // get the Lagrangian parameters
    LagrangianParams *p = (LagrangianParams *) params;
    double t0 = p->t0;
    doubleThreeVector pos0 = p->pos0;
    doubleThreeVector vel0 = p->vel0;

    // compute the Lagrangian
    std::array<doubleThreeVector, 2> newPosVel = CalculateNewPhaseSpaceCoords(t0, t, pos0, pos0, vel0, vel0);
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

    double epsabs = 1e-12;
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