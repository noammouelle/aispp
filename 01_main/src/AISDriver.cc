#include "AISDriver.hh"

AISDriver::AISDriver(AISParams params)
{
    // Set the frame rotation rate before anything else: the atom ensemble
    // constructors convert the requested frame velocities into the canonical
    // momenta the propagator works with, and that conversion needs Omega.
    bool rotatingPotential = (params.potentialType.rfind("rotating_", 0) == 0);
    bool rotationRequested = (params.rotationRate[0] != 0.0 ||
                              params.rotationRate[1] != 0.0 ||
                              params.rotationRate[2] != 0.0);
    if (rotationRequested && !rotatingPotential)
    {
        std::cerr << "Error: a non-zero 'rotation' was given but utype is '"
                  << params.potentialType << "', which ignores it. Use one of "
                  << "rotating_pot, rotating_linear_pot, rotating_quadratic_pot."
                  << std::endl;
        exit(1);
    }
    SetRotationRate(rotatingPotential ? params.rotationRate
                                      : doubleThreeVector{0.0, 0.0, 0.0});

    // create the atom ensemble
    if (params.usePhaseSpaceGrid) {
        atomEnsemble = std::make_unique<AISAtomEnsemble>(
            params.xGridMin,  params.xGridMax,  params.nGridX,
            params.yGridMin,  params.yGridMax,  params.nGridY,
            params.zGridMin,  params.zGridMax,  params.nGridZ,
            params.vxGridMin, params.vxGridMax, params.nGridVX,
            params.vyGridMin, params.vyGridMax, params.nGridVY,
            params.vzGridMin, params.vzGridMax, params.nGridVZ,
            params.initialState);
    } else {
        atomEnsemble = std::make_unique<AISAtomEnsemble>(params.nAtoms, params.cloudTransTemperature, params.cloudLongTemperature, params.cloudRadius,
                                                         params.initialPosition, params.initialVelocity,
                                                         params.seed, params.initialState);
    }

    // create the detector
    //detector = std::make_unique<AISDetector>(atomEnsemble, coherenceLength);
    std::unique_ptr<AISDetector> detector;

    // create the pointers to potential functions
    using potentialFunctionType = double(*)(const doubleThreeVector&, const doubleThreeVector&);
    using gradPotentialFunctionType = doubleThreeVector(*)(const doubleThreeVector&, const doubleThreeVector&);
    using hessianPotentialFunctionType = double3x3Matrix(*)(const doubleThreeVector&, const doubleThreeVector&);

    std::shared_ptr<potentialFunctionType> U;
    std::shared_ptr<gradPotentialFunctionType> dUdx, dUdp;
    std::shared_ptr<hessianPotentialFunctionType> d2Udxdx, d2Udxdp, d2Udpdp;

    // assign the potential function pointers
    if(params.potentialType == "zero_pot")
    {   U = std::make_shared<potentialFunctionType>(zeroU);
        dUdx = std::make_shared<gradPotentialFunctionType>(zeroGrad);
        dUdp = std::make_shared<gradPotentialFunctionType>(zeroGrad);
        d2Udxdx = std::make_shared<hessianPotentialFunctionType>(zeroHess);
        d2Udxdp = std::make_shared<hessianPotentialFunctionType>(zeroHess);
        d2Udpdp = std::make_shared<hessianPotentialFunctionType>(zeroHess);
    }
    else if(params.potentialType == "linear_pot")
    {
        U = std::make_shared<potentialFunctionType>(uniformGravityU);
        dUdx = std::make_shared<gradPotentialFunctionType>(uniformGravityGrad);
        dUdp = std::make_shared<gradPotentialFunctionType>(zeroGrad);
        d2Udxdx = std::make_shared<hessianPotentialFunctionType>(zeroHess);
        d2Udxdp = std::make_shared<hessianPotentialFunctionType>(zeroHess);
        d2Udpdp = std::make_shared<hessianPotentialFunctionType>(zeroHess);
    }
    else if(params.potentialType == "quadratic_pot")
    {
        U = std::make_shared<potentialFunctionType>(linearGravityU);
        dUdx = std::make_shared<gradPotentialFunctionType>(linearGravityGrad);
        dUdp = std::make_shared<gradPotentialFunctionType>(zeroGrad);
        d2Udxdx = std::make_shared<hessianPotentialFunctionType>(linearGravityHess);
        d2Udxdp = std::make_shared<hessianPotentialFunctionType>(zeroHess);
        d2Udpdp = std::make_shared<hessianPotentialFunctionType>(zeroHess);
    }
    // Rotating-frame variants. These differ from their inertial counterparts only
    // by the velocity-dependent Coriolis terms (dUdp and d2Udxdp); the
    // centrifugal force comes out of p^2/2m automatically. See AISPotentials.hh.
    else if(params.potentialType == "rotating_pot")
    {
        U = std::make_shared<potentialFunctionType>(rotatingU);
        dUdx = std::make_shared<gradPotentialFunctionType>(rotatingGrad);
        dUdp = std::make_shared<gradPotentialFunctionType>(rotationDUdp);
        d2Udxdx = std::make_shared<hessianPotentialFunctionType>(zeroHess);
        d2Udxdp = std::make_shared<hessianPotentialFunctionType>(rotationD2Udxdp);
        d2Udpdp = std::make_shared<hessianPotentialFunctionType>(zeroHess);
    }
    else if(params.potentialType == "rotating_linear_pot")
    {
        U = std::make_shared<potentialFunctionType>(rotatingUniformGravityU);
        dUdx = std::make_shared<gradPotentialFunctionType>(rotatingUniformGravityGrad);
        dUdp = std::make_shared<gradPotentialFunctionType>(rotationDUdp);
        d2Udxdx = std::make_shared<hessianPotentialFunctionType>(zeroHess);
        d2Udxdp = std::make_shared<hessianPotentialFunctionType>(rotationD2Udxdp);
        d2Udpdp = std::make_shared<hessianPotentialFunctionType>(zeroHess);
    }
    else if(params.potentialType == "rotating_quadratic_pot")
    {
        U = std::make_shared<potentialFunctionType>(rotatingLinearGravityU);
        dUdx = std::make_shared<gradPotentialFunctionType>(rotatingLinearGravityGrad);
        dUdp = std::make_shared<gradPotentialFunctionType>(rotationDUdp);
        d2Udxdx = std::make_shared<hessianPotentialFunctionType>(linearGravityHess);
        d2Udxdp = std::make_shared<hessianPotentialFunctionType>(rotationD2Udxdp);
        d2Udpdp = std::make_shared<hessianPotentialFunctionType>(zeroHess);
    }
    else
    {
        std::cerr << "Potential type "<<params.potentialType<<" not recognized. Exiting." << std::endl;
        exit(1);
    }

    // create the kinematic propagator
    kinematicPropagator = std::make_shared<AISKinematicPropagator>(U,dUdx,dUdp,d2Udxdx,d2Udxdp,d2Udpdp);
    kinematicPropagator->odeRelTol = params.gslKinRelError;
    kinematicPropagator->odeAbsTol = params.gslKinAbsError;
    kinematicPropagator->qagRelTol = params.gslQagRelError;
    kinematicPropagator->qagAbsTol = params.gslQagAbsError;
    kinematicPropagator->ultraFast = params.ultraFast;  

    __float128 omega_, omegaChirp_;
    doubleThreeVector k_, kChirp_;
    double rabiFreq_, phi0_, w0_, zLaser_, focalLength_, beamRadius_, baselineLength_;
    std::string beamType_;
    std::map<int, double> zernikeCoeffs_;

    for(int i = 0; i < params.rabiFrequencies.size(); ++i)
    {   
        if(params.wavefrontTypeVector[i] == "flat_square")
        {
            beamType_ = "flat_square";
        }
        else if(params.wavefrontTypeVector[i] == "gaussian")
        {
            beamType_ = "gaussian";
        }
        else if(params.wavefrontTypeVector[i] == "confocal")
        {
            beamType_ = "confocal";
        }
        else{
            std::cerr << "Wavefront type " << params.wavefrontTypeVector[i] << " not recognized. Exiting." << std::endl;
            exit(1);
        }
        // create the beam
        k_        = {params.kXVector[i], params.kYVector[i], params.kZVector[i]};
        omega_    = params.omegaVector[i];
        rabiFreq_ = params.rabiFrequencies[i];
        phi0_     = params.phi0[i];
        kChirp_   = {params.kXChirpVector[i], params.kYChirpVector[i], params.kZChirpVector[i]};
        omegaChirp_ = params.frequencyChirpVector[i];
        w0_ = params.waistVector[i];
        zLaser_ = params.zLaserVector[i];
        focalLength_ = params.focalLengthVector[i];
        beamRadius_ = params.beamRadiusVector[i];
        baselineLength_ = params.baselineLengthVector[i];


        // get the zernike coeffs
        for(const auto& [key, value] : params.zernikeCoeff)
        {
            zernikeCoeffs_[key] = value[i];
        }

        std::shared_ptr<AISLaserBeam> beam = std::make_shared<AISLaserBeam>(k_,omega_,rabiFreq_, phi0_);

        beam->SetKChirp(kChirp_);
        beam->SetFrequencyChirp(omegaChirp_);
        beam->SetW0(w0_);
        beam->SetZLaser(zLaser_);
        beam->SetFocalLength(focalLength_);        
        beam->SetZernikeCoeffs(zernikeCoeffs_);
        beam->SetBeamRadius(beamRadius_);
        beam->SetBeamType(beamType_);
        beam->SetBaselineLength(baselineLength_);

        // create the propagator
        auto pulsePropagator = std::make_shared<AISPulsePropagator>(beam, params.initialPulseTimes[i], params.finalPulseTimes[i], kinematicPropagator,
                                                                   params.amplitudeThreshold);
        pulsePropagator->SetIgnoreDetuning(params.ignoreDetuning);
        pulsePropagator->relTol = params.gslPulseRelError;
        pulsePropagator->absTol = params.gslPulseAbsError;
        pulsePropagator->ultraFast = params.ultraFast;

        pulsePropagator->SetUsePathSelection(params.usePathSelection);
        pulsePropagator->SetPathsToSimulate(params.pathsToSimulate);

        pulsePropagator->SetUseStaticApprox(params.useStaticApprox);

        // add the propagator to the list
        pulsePropagators.push_back(pulsePropagator);
    }

    // store the params
    this->params = params;
}

AISDriver::~AISDriver()
{
}

void AISDriver::Run()
{
    // 1. If required, run the path finder sequence
    if(params.useMcBranching or params.useDetVolSelection)
    {
        std::cout<<"Running path finder sequence..."<<std::endl;
        // start the timer
        auto start = std::chrono::high_resolution_clock::now();
        RunPathFinder();
        // stop the timer
        auto stop = std::chrono::high_resolution_clock::now();

        // set the interferingPaths variable in the pulse propagators
        for(int i = 0; i < pulsePropagators.size(); ++i)
        {
            if(params.useMcBranching)
            {
            this->pulsePropagators[i]->SetInterferingPaths(interferingPaths);
            this->pulsePropagators[i]->SetUseMcBranching(true);
            }
            if(params.useDetVolSelection)
            {
                this->pulsePropagators[i]->SetUseDetVolSelection(true);
            }
        }
    }
    // if using detection volume selection, set the detectable paths for each wavepacket
    if(params.useDetVolSelection)
    {
        for(int i = 0; i < params.nAtoms; ++i)
        {
            std::unique_ptr<AISWavePacket>& currentWavePacket = this->atomEnsemble->GetAtom(i)->GetWavePacket(0);
            // compare frame velocities: params.initialVelocity is a frame velocity
            doubleThreeVector wpPos = currentWavePacket->GetPosition();
            currentWavePacket->SetDetectablePaths(
                GetDetectablePaths(wpPos, frameFromCanonicalVelocity(wpPos, currentWavePacket->GetVelocity())));
        }
    }

    // 2. Run the main sequence
    std::cout<<"Running main sequence..."<<std::endl;
    RunAll();

    // printout the total number of wavepackets
    int n = 0;
    for(int i = 0; i < atomEnsemble->GetNumberOfAtoms(); ++i)
    {
        n += atomEnsemble->GetAtom(i)->GetNumberOfWavePackets();
    }
    std::cout<<"Total number of wavepackets: "<<n<<std::endl;
}

void AISDriver::RunAll()
{
    bool traj = params.printTrajectory;
    if (traj && atomEnsemble->GetNumberOfAtoms() > 10)
        std::cerr << "Warning: printtrajectory=1 with >10 atoms — "
                     "trajectory files may be large." << std::endl;

    if (traj) RecordSnapshot(0.0, "initial");

    // initial propagation, if the initial pulse time is not zero
    if(params.initialPulseTimes[0] != 0.0q)
    {
        kinematicPropagator->SetAddEnergyPhase(false);
        kinematicPropagator->PropagateEnsemble(atomEnsemble, params.initialPulseTimes[0]);
        if (traj) RecordSnapshot(static_cast<double>(params.initialPulseTimes[0]),
                                  "pulse_0_start");
    }

    // loop over the pulse propagators, and kinematic propagators in between
    for(int i = 0; i < pulsePropagators.size(); ++i)
    {
        // pulse propagation
        pulsePropagators[i]->PropagateEnsemble(atomEnsemble);
        if (traj) RecordSnapshot(static_cast<double>(params.finalPulseTimes[i]),
                                  "pulse_" + std::to_string(i) + "_end");

        // kinematic propagation (unless it is the last pulse)
        if(i < pulsePropagators.size() - 1)
        {
            kinematicPropagator->SetAddEnergyPhase(true);
            kinematicPropagator->PropagateEnsemble(atomEnsemble, params.initialPulseTimes[i + 1]);
            if (traj) RecordSnapshot(static_cast<double>(params.initialPulseTimes[i+1]),
                                      "pulse_" + std::to_string(i+1) + "_start");
        }
    }

    // final propagation
    if(params.detectionTime != params.finalPulseTimes.back())
    {
        kinematicPropagator->SetAddEnergyPhase(false);
        kinematicPropagator->PropagateEnsemble(atomEnsemble, params.detectionTime);
        if (traj) RecordSnapshot(static_cast<double>(params.detectionTime), "detection");
    }
}

void AISDriver::RunPathFinder()
{
    AISParams pathFinderParams = params;
    // set position and velocity spread to 0
    pathFinderParams.cloudRadius = 0.0;
    pathFinderParams.cloudTransTemperature = 0.0;
    pathFinderParams.cloudLongTemperature = 0.0;
    // set the number of atoms to 1
    pathFinderParams.nAtoms = 1;
    // set the useMcBranching flag to false
    pathFinderParams.useMcBranching = false;
    // use large cutoff, as we are only interested in the central (interfering) wavepackets
    if(pathFinderParams.amplitudeThreshold < 0.001)
    {
        pathFinderParams.amplitudeThreshold = 0.001; // keep the same as in the main sequence otherwise
    }

    // create another driver
    std::unique_ptr<AISDriver> pathFinder = std::make_unique<AISDriver>(pathFinderParams);

    // run the path finder
    pathFinder->RunAll();

    // create the detector
    detector = std::make_unique<AISDetector>(pathFinder->atomEnsemble, pathFinderParams.coherenceLength);

    // get the paths of the interfering atoms
    interferingPaths = detector->GetPortFrame(0)->GetInterferingPaths();

    // get the final positions of all wavepackets
    for(int i = 0; i < pathFinder->atomEnsemble->GetAtom(0)->GetNumberOfWavePackets(); ++i)
    {
        std::unique_ptr<AISWavePacket>& currentWavePacket = pathFinder->atomEnsemble->GetAtom(0)->GetWavePacket(i);
        pathToFinalPosMap[currentWavePacket->GetPath()] = currentWavePacket->GetPosition();
    }
}

void AISDriver::Detect()
{   
    // init the detector
    detector = std::make_unique<AISDetector>(atomEnsemble, params.coherenceLength, params.ultraFast);
    sampledPortIndices = detector->SampleAllPorts();
}

void AISDriver::WriteDetectedAtomsToFile(std::string filename)
{
    // open a h5 file
    H5::H5File file(filename, H5F_ACC_TRUNC);
    int numSamples = detector->GetNumberOfSamples();

    // create the data vectors
    std::vector<doubleThreeVector> positions(numSamples);
    std::vector<doubleThreeVector> velocities(numSamples);
    std::vector<int> states(numSamples);
    std::vector<int> interferingFlag(numSamples);
    std::vector<double> phaseShifts(numSamples);
    std::vector<double> phaseShiftErrors(numSamples);

    // loop over the samples to fill in the data
    int sampleIndex = 0;
    for(int portFrameIndex = 0; portFrameIndex < detector->GetNumberOfPortFrames(); portFrameIndex++)
    {    
        int sampledPortIndex = sampledPortIndices[portFrameIndex];
        if(sampledPortIndex != -1)
        {
            // get the port frame
            std::unique_ptr<AISPortFrame>& portFrame = detector->GetPortFrame(portFrameIndex);
            // get the port
            std::unique_ptr<AISPort>& port = portFrame->GetPort(sampledPortIndex);
            // get the data
            doubleThreeVector currentPosition = port->position;
            // report the frame velocity, not the canonical p/m
            doubleThreeVector currentVelocity = frameFromCanonicalVelocity(currentPosition,
                                                                           port->velocity);
            int currentState = port->state;
            int currentInterferingFlag = static_cast<int>(port->interfering);
            double currentPhaseShift = port->phaseShift;
            double currentPhaseShiftError = port->phaseShiftError;
            // fill in the data
            for(int i = 0; i < 3; ++i)
            {
                positions[sampleIndex][i] = currentPosition[i];
                velocities[sampleIndex][i] = currentVelocity[i];
            }
            states[sampleIndex] = currentState;
            interferingFlag[sampleIndex] = currentInterferingFlag;
            phaseShifts[sampleIndex] = currentPhaseShift;
            phaseShiftErrors[sampleIndex] = currentPhaseShiftError;

            sampleIndex++;
        }        
    }
    // define the dimensionality of the datasets
    hsize_t numSamplesHsize = numSamples;
    hsize_t dim_positions[2] = {numSamplesHsize, 3};
    hsize_t dim_velocities[2] = {numSamplesHsize, 3};
    hsize_t dim_states[1] = {numSamplesHsize};
    hsize_t dim_interferingFlag[1] = {numSamplesHsize};
    hsize_t dim_phaseShifts[1] = {numSamplesHsize};
    hsize_t dim_phaseShiftErrors[1] = {numSamplesHsize};
    // create the dataspace
    H5::DataSpace dataspace_positions(2, dim_positions);
    H5::DataSpace dataspace_velocities(2, dim_velocities);
    H5::DataSpace dataspace_states(1, dim_states);
    H5::DataSpace dataspace_interferingFlag(1, dim_interferingFlag);
    H5::DataSpace dataspace_phaseShifts(1, dim_phaseShifts);
    H5::DataSpace dataspace_phaseShiftErrors(1, dim_phaseShiftErrors);
    // create the datasets
    H5::DataSet dataset_positions = file.createDataSet("positions", H5::PredType::NATIVE_DOUBLE, dataspace_positions);
    H5::DataSet dataset_velocities = file.createDataSet("velocities", H5::PredType::NATIVE_DOUBLE, dataspace_velocities);
    H5::DataSet dataset_states = file.createDataSet("states", H5::PredType::NATIVE_INT, dataspace_states);
    H5::DataSet dataset_interferingFlag = file.createDataSet("interferingFlag", H5::PredType::NATIVE_INT, dataspace_interferingFlag);
    H5::DataSet dataset_phaseShifts = file.createDataSet("phaseShifts", H5::PredType::NATIVE_DOUBLE, dataspace_phaseShifts);
    H5::DataSet dataset_phaseShiftErrors = file.createDataSet("phaseShiftErrors", H5::PredType::NATIVE_DOUBLE, dataspace_phaseShiftErrors);
    // write the data
    dataset_positions.write(positions.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_velocities.write(velocities.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_states.write(states.data(), H5::PredType::NATIVE_INT);
    dataset_interferingFlag.write(interferingFlag.data(), H5::PredType::NATIVE_INT);
    dataset_phaseShifts.write(phaseShifts.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_phaseShiftErrors.write(phaseShiftErrors.data(), H5::PredType::NATIVE_DOUBLE);

    WriteRotationRate(file);

    // close the file
    file.close();
}

void AISDriver::WriteWavePacketsToFile(std::string fName)
{
    writeAtomEnsembleToFile(fName, atomEnsemble);
}

void AISDriver::WritePortsToFile(std::string fName)
{
    // write the state, position, velocity, interference flag and probability of each port, regardless of atom

    // open a h5 file
    H5::H5File file(fName, H5F_ACC_TRUNC);

    // count number of ports
    int numPorts=0;
    for(int i=0; i < detector->GetNumberOfPortFrames(); i++)
    {
        std::unique_ptr<AISPortFrame>& portFrame = detector->GetPortFrame(i);
        for(int j=0; j < portFrame->GetNumberOfPorts(); j++)
        {
            numPorts++;
        }
    }

    // create the data vectors
    std::vector<doubleThreeVector> positions(numPorts);
    std::vector<doubleThreeVector> velocities(numPorts);
    std::vector<int> states(numPorts);
    std::vector<int> interferingFlag(numPorts);
    std::vector<double> probabilities(numPorts);
    std::vector<double> phaseShifts(numPorts);
    std::vector<double> phaseShiftErrors(numPorts);

    int index = 0;

    // loop over the ports to fill in the data
    for(int portFrameIndex = 0; portFrameIndex < detector->GetNumberOfPortFrames(); portFrameIndex++)
    {    
        std::unique_ptr<AISPortFrame>& portFrame = detector->GetPortFrame(portFrameIndex);
        for(int portIndex = 0; portIndex < portFrame->GetNumberOfPorts(); portIndex++)
        {
            // get the port
            std::unique_ptr<AISPort>& port = portFrame->GetPort(portIndex);
            // get the data
            doubleThreeVector currentPosition = port->position;
            // report the frame velocity, not the canonical p/m
            doubleThreeVector currentVelocity = frameFromCanonicalVelocity(currentPosition,
                                                                           port->velocity);
            int currentState = port->state;
            int currentInterferingFlag = static_cast<int>(port->interfering);
            double currentProbability = abs(port->probabilityAmplitude);
            double currentPhaseShift = port->phaseShift;
            double currentPhaseShiftError = port->phaseShiftError;

            // fill in the data
            for(int i = 0; i < 3; ++i)
            {
                positions[index][i] = currentPosition[i];
                velocities[index][i] = currentVelocity[i];
            }
            states[index] = currentState;
            interferingFlag[index] = currentInterferingFlag;
            probabilities[index] = currentProbability;
            phaseShifts[index] = currentPhaseShift;
            phaseShiftErrors[index] = currentPhaseShiftError;

            index++;
        }        
    }
    // define the dimensionality of the datasets
    hsize_t numSPortsHsize = numPorts;
    hsize_t dim_positions[2] = {numSPortsHsize, 3};
    hsize_t dim_velocities[2] = {numSPortsHsize, 3};
    hsize_t dim_states[1] = {numSPortsHsize};
    hsize_t dim_interferingFlag[1] = {numSPortsHsize};
    hsize_t dim_probabilities[1] = {numSPortsHsize};
    hsize_t dim_phaseShifts[1] = {numSPortsHsize};
    hsize_t dim_phaseShiftErrors[1] = {numSPortsHsize};
    // create the dataspace
    H5::DataSpace dataspace_positions(2, dim_positions);
    H5::DataSpace dataspace_velocities(2, dim_velocities);
    H5::DataSpace dataspace_states(1, dim_states);
    H5::DataSpace dataspace_interferingFlag(1, dim_interferingFlag);
    H5::DataSpace dataspace_probabilities(1, dim_probabilities);
    H5::DataSpace dataspace_phaseShifts(1, dim_phaseShifts);
    H5::DataSpace dataspace_phaseShiftErrors(1, dim_phaseShiftErrors);
    // create the datasets
    H5::DataSet dataset_positions = file.createDataSet("positions", H5::PredType::NATIVE_DOUBLE, dataspace_positions);
    H5::DataSet dataset_velocities = file.createDataSet("velocities", H5::PredType::NATIVE_DOUBLE, dataspace_velocities);
    H5::DataSet dataset_states = file.createDataSet("states", H5::PredType::NATIVE_INT, dataspace_states);
    H5::DataSet dataset_interferingFlag = file.createDataSet("interferingFlag", H5::PredType::NATIVE_INT, dataspace_interferingFlag);
    H5::DataSet dataset_probabilities = file.createDataSet("probabilities", H5::PredType::NATIVE_DOUBLE, dataspace_probabilities);
    H5::DataSet dataset_phaseShifts = file.createDataSet("phaseShifts", H5::PredType::NATIVE_DOUBLE, dataspace_phaseShifts);
    H5::DataSet dataset_phaseShiftErrors = file.createDataSet("phaseShiftErrors", H5::PredType::NATIVE_DOUBLE, dataspace_phaseShiftErrors);
    // write the data
    dataset_positions.write(positions.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_velocities.write(velocities.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_states.write(states.data(), H5::PredType::NATIVE_INT);
    dataset_interferingFlag.write(interferingFlag.data(), H5::PredType::NATIVE_INT);
    dataset_probabilities.write(probabilities.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_phaseShifts.write(phaseShifts.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_phaseShiftErrors.write(phaseShiftErrors.data(), H5::PredType::NATIVE_DOUBLE);

    WriteRotationRate(file);

    // close the file
    file.close();
}

void AISDriver::WritePhaseSpaceMapToFile(std::string fName)
{
    // detector must already be populated by Detect() before calling this.
    int nAtoms = atomEnsemble->GetNumberOfAtoms();

    // Count total ports across all atoms
    int totalPorts = 0;
    for (int ai = 0; ai < nAtoms; ++ai)
        totalPorts += detector->GetPortFrame(ai)->GetNumberOfPorts();

    // One row per port
    std::vector<double> initPosFlat(totalPorts * 3), initVelFlat(totalPorts * 3);
    std::vector<double> finalPosFlat(totalPorts * 3), finalVelFlat(totalPorts * 3);
    std::vector<double> phaseShifts(totalPorts, 0.0), phaseShiftErrors(totalPorts, 0.0);
    std::vector<int>    states(totalPorts), interferingFlag(totalPorts), atomIndices(totalPorts);
    std::vector<double> amp0(totalPorts, 0.0), amp1(totalPorts, 0.0);
    std::vector<std::string> path0Strs(totalPorts, ""), path1Strs(totalPorts, "");

    int idx = 0;
    for (int ai = 0; ai < nAtoms; ++ai)
    {
        auto& atom      = atomEnsemble->GetAtom(ai);
        auto& portFrame = detector->GetPortFrame(ai);

        // Initial phase-space coords — same for every wavepacket of this atom
        doubleThreeVector p0 = atom->GetWavePacket(0)->GetPos0();
        // report frame velocities, not canonical p/m
        doubleThreeVector v0 = frameFromCanonicalVelocity(p0, atom->GetWavePacket(0)->GetVel0());

        for (int pi = 0; pi < portFrame->GetNumberOfPorts(); ++pi)
        {
            auto& port = portFrame->GetPort(pi);

            for (int k = 0; k < 3; ++k)
            {
                initPosFlat [idx*3+k] = p0[k];
                initVelFlat [idx*3+k] = v0[k];
                finalPosFlat[idx*3+k] = port->position[k];
                finalVelFlat[idx*3+k] = frameFromCanonicalVelocity(port->position,
                                                                    port->velocity)[k];
            }

            phaseShifts     [idx] = port->phaseShift;
            phaseShiftErrors[idx] = port->phaseShiftError;
            states          [idx] = port->state;
            interferingFlag [idx] = static_cast<int>(port->interfering);
            atomIndices     [idx] = ai;

            // WP 0 — always present
            auto& wp0 = atom->GetWavePacket(port->wavePacketIndices[0]);
            amp0[idx]      = wp0->GetAmplitude();
            path0Strs[idx] = wp0->GetPath();

            // WP 1 — only for interfering (2-wavepacket) ports
            if (port->interfering && port->getNumberOfWavePackets() >= 2)
            {
                auto& wp1  = atom->GetWavePacket(port->wavePacketIndices[1]);
                amp1[idx]      = wp1->GetAmplitude();
                path1Strs[idx] = wp1->GetPath();
            }

            ++idx;
        }
    }

    H5::H5File file(fName, H5F_ACC_TRUNC);

    hsize_t n = static_cast<hsize_t>(totalPorts);
    hsize_t dim1D[1] = {n};
    hsize_t dim2D[2] = {n, 3};

    auto write1D_double = [&](const std::string& name, const std::vector<double>& data) {
        H5::DataSpace ds(1, dim1D);
        file.createDataSet(name, H5::PredType::NATIVE_DOUBLE, ds).write(data.data(), H5::PredType::NATIVE_DOUBLE);
    };
    auto write1D_int = [&](const std::string& name, const std::vector<int>& data) {
        H5::DataSpace ds(1, dim1D);
        file.createDataSet(name, H5::PredType::NATIVE_INT, ds).write(data.data(), H5::PredType::NATIVE_INT);
    };
    auto write2D_double = [&](const std::string& name, const std::vector<double>& data) {
        H5::DataSpace ds(2, dim2D);
        file.createDataSet(name, H5::PredType::NATIVE_DOUBLE, ds).write(data.data(), H5::PredType::NATIVE_DOUBLE);
    };
    auto write1D_str = [&](const std::string& name, const std::vector<std::string>& strs) {
        H5::StrType strType(H5::PredType::C_S1, H5T_VARIABLE);
        H5::DataSpace ds(1, dim1D);
        std::vector<const char*> ptrs(strs.size());
        for (size_t i = 0; i < strs.size(); ++i) ptrs[i] = strs[i].c_str();
        file.createDataSet(name, strType, ds).write(ptrs.data(), strType);
    };

    write2D_double("initial_positions",  initPosFlat);
    write2D_double("initial_velocities", initVelFlat);
    write2D_double("final_positions",    finalPosFlat);
    write2D_double("final_velocities",   finalVelFlat);
    write1D_double("phase_shifts",       phaseShifts);
    write1D_double("phase_shift_errors", phaseShiftErrors);
    write1D_int   ("states",             states);
    write1D_int   ("is_interfering",     interferingFlag);
    write1D_int   ("atom_indices",       atomIndices);
    write1D_double("amp0",               amp0);
    write1D_double("amp1",               amp1);
    write1D_str   ("path0",             path0Strs);
    write1D_str   ("path1",             path1Strs);

    WriteRotationRate(file);

    file.close();
}

void AISDriver::WriteRotationRate(H5::H5File& file)
{
    hsize_t dim[1] = {3};
    H5::DataSpace ds(1, dim);
    // use the rate the physics actually ran with, not the requested one
    std::array<double,3> rot = {gRotationRate[0], gRotationRate[1], gRotationRate[2]};
    file.createDataSet("rotation", H5::PredType::NATIVE_DOUBLE, ds)
        .write(rot.data(), H5::PredType::NATIVE_DOUBLE);
}

void AISDriver::RecordSnapshot(double time, const std::string& label)
{
    int snapIdx = static_cast<int>(fTrajTimes.size());
    fTrajTimes.push_back(time);
    fTrajLabels.push_back(label);

    int nAtoms = atomEnsemble->GetNumberOfAtoms();
    for (int ai = 0; ai < nAtoms; ++ai)
    {
        auto& atom = atomEnsemble->GetAtom(ai);
        for (int wi = 0; wi < atom->GetNumberOfWavePackets(); ++wi)
        {
            auto& wp = atom->GetWavePacket(wi);
            doubleThreeVector pos = wp->GetPosition();
            // snapshots record the frame velocity so that the trajectory traces
            // aispy reconstructs are the physical ones
            doubleThreeVector vel = frameFromCanonicalVelocity(pos, wp->GetVelocity());
            fTrajSnapIdx.push_back(snapIdx);
            fTrajAtomIdx.push_back(ai);
            fTrajPaths.push_back(wp->GetPath());
            fTrajStates.push_back(wp->GetState());
            fTrajAmplitudes.push_back(wp->GetAmplitude());
            fTrajPositions.push_back({pos[0], pos[1], pos[2]});
            fTrajVelocities.push_back({vel[0], vel[1], vel[2]});
        }
    }
}

void AISDriver::WriteTrajectoryToFile(std::string fName)
{
    int nSnap = static_cast<int>(fTrajTimes.size());
    int nRec  = static_cast<int>(fTrajSnapIdx.size());

    if (nRec == 0)
    {
        std::cerr << "WriteTrajectoryToFile: no snapshots recorded. "
                     "Did you set printtrajectory 1?" << std::endl;
        return;
    }

    H5::H5File file(fName, H5F_ACC_TRUNC);

    // ── snapshot-level datasets ───────────────────────────────────────────────
    {
        hsize_t dim[1] = {static_cast<hsize_t>(nSnap)};
        H5::DataSpace ds(1, dim);
        file.createDataSet("snapshot_times", H5::PredType::NATIVE_DOUBLE, ds)
            .write(fTrajTimes.data(), H5::PredType::NATIVE_DOUBLE);

        H5::StrType strType(H5::PredType::C_S1, H5T_VARIABLE);
        std::vector<const char*> labelPtrs(nSnap);
        for (int i = 0; i < nSnap; ++i) labelPtrs[i] = fTrajLabels[i].c_str();
        file.createDataSet("snapshot_labels", strType, ds)
            .write(labelPtrs.data(), strType);
    }

    // ── per-record datasets ───────────────────────────────────────────────────
    hsize_t n   = static_cast<hsize_t>(nRec);
    hsize_t n3[2] = {n, 3};
    H5::DataSpace ds1D(1, &n);
    H5::DataSpace ds2D(2, n3);

    auto write1D_int = [&](const std::string& name, const std::vector<int>& v) {
        file.createDataSet(name, H5::PredType::NATIVE_INT, ds1D)
            .write(v.data(), H5::PredType::NATIVE_INT);
    };
    auto write1D_dbl = [&](const std::string& name, const std::vector<double>& v) {
        file.createDataSet(name, H5::PredType::NATIVE_DOUBLE, ds1D)
            .write(v.data(), H5::PredType::NATIVE_DOUBLE);
    };
    auto write1D_str = [&](const std::string& name, const std::vector<std::string>& v) {
        H5::StrType strType(H5::PredType::C_S1, H5T_VARIABLE);
        std::vector<const char*> ptrs(v.size());
        for (size_t i = 0; i < v.size(); ++i) ptrs[i] = v[i].c_str();
        file.createDataSet(name, strType, ds1D).write(ptrs.data(), strType);
    };

    write1D_int("snapshot_idx", fTrajSnapIdx);
    write1D_int("atom_indices", fTrajAtomIdx);
    write1D_str("paths",        fTrajPaths);
    write1D_int("states",       fTrajStates);
    write1D_dbl("amplitudes",   fTrajAmplitudes);

    // Flatten 3-vectors
    std::vector<double> posFlat(nRec*3), velFlat(nRec*3);
    for (int i = 0; i < nRec; ++i)
    {
        for (int k = 0; k < 3; ++k)
        {
            posFlat[i*3+k] = fTrajPositions[i][k];
            velFlat[i*3+k] = fTrajVelocities[i][k];
        }
    }
    file.createDataSet("positions",  H5::PredType::NATIVE_DOUBLE, ds2D)
        .write(posFlat.data(), H5::PredType::NATIVE_DOUBLE);
    file.createDataSet("velocities", H5::PredType::NATIVE_DOUBLE, ds2D)
        .write(velFlat.data(), H5::PredType::NATIVE_DOUBLE);

    // Frame rotation rate, so that aispy can reconstruct the free flight between
    // snapshots without having to re-read the input file.
    WriteRotationRate(file);

    file.close();
    std::cout << "Trajectory: " << nSnap << " snapshots, "
              << nRec << " records → " << fName << std::endl;
}

std::vector<std::string> AISDriver::GetDetectablePaths(doubleThreeVector pos0, doubleThreeVector vel0)
{
    std::vector<std::string> detectablePaths;
    // compute the delta with central position/velocity of the cloud
    doubleThreeVector cloudPos0 = params.initialPosition;
    doubleThreeVector cloudVel0 = params.initialVelocity;
    doubleThreeVector dPos = {pos0[0] - cloudPos0[0], pos0[1] - cloudPos0[1], pos0[2] - cloudPos0[2]};
    doubleThreeVector dVel = {vel0[0] - cloudVel0[0], vel0[1] - cloudVel0[1], vel0[2] - cloudVel0[2]};
    // loop over the path map
    for(auto const& pathToFinalPos : pathToFinalPosMap)
    {
        // get the path
        std::string path = pathToFinalPos.first;
        // get the final position
        doubleThreeVector finalPos = pathToFinalPos.second;
        
        // linearly propagate the final position of the wavepacket to the detection time
        doubleThreeVector translatedFinalPos = finalPos;
        for(int i = 0; i < 3; ++i)
        {
            translatedFinalPos[i] += dPos[i] + dVel[i] * params.detectionTime;
        }

        // now check if the translated final position is within the detection volume
        if( translatedFinalPos[0] > params.xDetMin && translatedFinalPos[0] < params.xDetMax &&
            translatedFinalPos[1] > params.yDetMin && translatedFinalPos[1] < params.yDetMax &&
            translatedFinalPos[2] > params.zDetMin && translatedFinalPos[2] < params.zDetMax)
        {
            detectablePaths.push_back(path);
        }
    }

    return detectablePaths;
}