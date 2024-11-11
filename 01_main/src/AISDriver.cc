#include "AISDriver.hh"

AISDriver::AISDriver(AISParams params)
{   
    // create the atom ensemble
    atomEnsemble = std::make_unique<AISAtomEnsemble>(params.nAtoms, params.cloudTemperature, params.cloudRadius,
                                                     params.initialPosition, params.initialVelocity);
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
    else
    {
        std::cerr << "Potential type "<<params.potentialType<<" not recognized. Exiting." << std::endl;
        exit(1);
    }

    // create the kinematic propagator
    kinematicPropagator = std::make_shared<AISKinematicPropagator>(U,dUdx,dUdp,d2Udxdx,d2Udxdp,d2Udpdp);

    // create the pulse propagators
    using wavefrontFunctionType = double(*)(const doubleThreeVector&);
    using delWavefrontFunctionType = doubleThreeVector(*)(const doubleThreeVector&);
    using rabifreqFunctionType  = double(*)(const doubleThreeVector&, const __float128&,const  __float128&);

    std::shared_ptr<wavefrontFunctionType> wff;
    std::shared_ptr<delWavefrontFunctionType> dwff;
    std::shared_ptr<rabifreqFunctionType> rff;

    __float128 omega_;
    doubleThreeVector k_;
    double rabiFreq_, phi0_;

    for(int i = 0; i < params.rabiFrequencies.size(); ++i)
    {   
        if(params.wavefrontTypeVector[i] == "flat_square")
        {
            wff = std::make_shared<wavefrontFunctionType>(flatWavefront);
            dwff = std::make_shared<delWavefrontFunctionType>(flatGradientWavefront);
            rff = std::make_shared<rabifreqFunctionType>(flatSquareEnvelope);
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

        std::shared_ptr<AISLaserBeam> beam = std::make_shared<AISLaserBeam>(k_,omega_,rabiFreq_, phi0_,
                                                                            wff,dwff,rff);

        // create the propagator
        auto pulsePropagator = std::make_shared<AISPulsePropagator>(beam, params.initialPulseTimes[i], params.finalPulseTimes[i], kinematicPropagator,
                                                                   params.amplitudeThreshold);

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
    __float128 t0, t1;
    // initial propagation, if the initial pulse time is not zero
    if(params.initialPulseTimes[0] != 0.0q)
    {
        kinematicPropagator->SetAddEnergyPhase(false);
        kinematicPropagator->PropagateEnsemble(atomEnsemble, params.initialPulseTimes[0]);
    }

    // loop over the pulse propagators, and kinematic propagators in between
    for(int i = 0; i < pulsePropagators.size(); ++i)
    {
        // pulse propagation
        pulsePropagators[i]->PropagateEnsemble(atomEnsemble);

        // kinematic propagation (unless it is the last pulse)
        if(i < pulsePropagators.size() - 1)
        {
            kinematicPropagator->SetAddEnergyPhase(true);
            kinematicPropagator->PropagateEnsemble(atomEnsemble, params.initialPulseTimes[i + 1]); // propagate to the start of the next pulse
        }
    }

    // final propagation
    if(params.detectionTime != params.finalPulseTimes.back())
    {
        kinematicPropagator->SetAddEnergyPhase(false);
        kinematicPropagator->PropagateEnsemble(atomEnsemble, params.detectionTime);
    }
}

void AISDriver::Detect()
{   
    // three times the De Broglie wavelength of the cloud
    double coherenceLength = 0.5 * hbar / sqrt(massSr87 * kB * params.cloudTemperature);
    // init the detector
    detector = std::make_unique<AISDetector>(atomEnsemble, 10 * coherenceLength); // 10 times the coherence length is the interference limit
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
            doubleThreeVector currentVelocity = port->velocity;
            int currentState = port->state;
            int currentInterferingFlag = static_cast<int>(port->interfering);
            double currentPhaseShift = port->phaseShift;
            // fill in the data
            for(int i = 0; i < 3; ++i)
            {
                positions[sampleIndex][i] = currentPosition[i];
                velocities[sampleIndex][i] = currentVelocity[i];
            }
            states[sampleIndex] = currentState;
            interferingFlag[sampleIndex] = currentInterferingFlag;
            phaseShifts[sampleIndex] = currentPhaseShift;

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
    // create the dataspace
    H5::DataSpace dataspace_positions(2, dim_positions);
    H5::DataSpace dataspace_velocities(2, dim_velocities);
    H5::DataSpace dataspace_states(1, dim_states);
    H5::DataSpace dataspace_interferingFlag(1, dim_interferingFlag);
    H5::DataSpace dataspace_phaseShifts(1, dim_phaseShifts);
    // create the datasets
    H5::DataSet dataset_positions = file.createDataSet("positions", H5::PredType::NATIVE_DOUBLE, dataspace_positions);
    H5::DataSet dataset_velocities = file.createDataSet("velocities", H5::PredType::NATIVE_DOUBLE, dataspace_velocities);
    H5::DataSet dataset_states = file.createDataSet("states", H5::PredType::NATIVE_INT, dataspace_states);
    H5::DataSet dataset_interferingFlag = file.createDataSet("interferingFlag", H5::PredType::NATIVE_INT, dataspace_interferingFlag);
    H5::DataSet dataset_phaseShifts = file.createDataSet("phaseShifts", H5::PredType::NATIVE_DOUBLE, dataspace_phaseShifts);
    // write the data
    dataset_positions.write(positions.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_velocities.write(velocities.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_states.write(states.data(), H5::PredType::NATIVE_INT);
    dataset_interferingFlag.write(interferingFlag.data(), H5::PredType::NATIVE_INT);
    dataset_phaseShifts.write(phaseShifts.data(), H5::PredType::NATIVE_DOUBLE);

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
            doubleThreeVector currentVelocity = port->velocity;
            int currentState = port->state;
            int currentInterferingFlag = static_cast<int>(port->interfering);
            double currentProbability = abs(port->probabilityAmplitude) * abs(port->probabilityAmplitude);
            // fill in the data
            for(int i = 0; i < 3; ++i)
            {
                positions[portIndex][i] = currentPosition[i];
                velocities[portIndex][i] = currentVelocity[i];
            }
            states[portIndex] = currentState;
            interferingFlag[portIndex] = currentInterferingFlag;
            probabilities[portIndex] = currentProbability;
        }        
    }
    // define the dimensionality of the datasets
    hsize_t numSPortsHsize = numPorts;
    hsize_t dim_positions[2] = {numSPortsHsize, 3};
    hsize_t dim_velocities[2] = {numSPortsHsize, 3};
    hsize_t dim_states[1] = {numSPortsHsize};
    hsize_t dim_interferingFlag[1] = {numSPortsHsize};
    hsize_t dim_probabilities[1] = {numSPortsHsize};
    // create the dataspace
    H5::DataSpace dataspace_positions(2, dim_positions);
    H5::DataSpace dataspace_velocities(2, dim_velocities);
    H5::DataSpace dataspace_states(1, dim_states);
    H5::DataSpace dataspace_interferingFlag(1, dim_interferingFlag);
    H5::DataSpace dataspace_probabilities(1, dim_probabilities);
    // create the datasets
    H5::DataSet dataset_positions = file.createDataSet("positions", H5::PredType::NATIVE_DOUBLE, dataspace_positions);
    H5::DataSet dataset_velocities = file.createDataSet("velocities", H5::PredType::NATIVE_DOUBLE, dataspace_velocities);
    H5::DataSet dataset_states = file.createDataSet("states", H5::PredType::NATIVE_INT, dataspace_states);
    H5::DataSet dataset_interferingFlag = file.createDataSet("interferingFlag", H5::PredType::NATIVE_INT, dataspace_interferingFlag);
    H5::DataSet dataset_probabilities = file.createDataSet("probabilities", H5::PredType::NATIVE_DOUBLE, dataspace_probabilities);
    // write the data
    dataset_positions.write(positions.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_velocities.write(velocities.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_states.write(states.data(), H5::PredType::NATIVE_INT);
    dataset_interferingFlag.write(interferingFlag.data(), H5::PredType::NATIVE_INT);
    dataset_probabilities.write(probabilities.data(), H5::PredType::NATIVE_DOUBLE);

    // close the file
    file.close();
}
