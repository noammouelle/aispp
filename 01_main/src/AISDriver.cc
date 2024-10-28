#include "AISDriver.hh"

AISDriver::AISDriver(AISParams params)
{   
    // three times the De Broglie wavelength of the cloud
    double coherenceLength = 3 * sqrt(2 * pi) * hbar / sqrt(massSr87 * kB * params.cloudTemperature);
    // create the atom ensemble
    atomEnsemble = new AISAtomEnsemble(params.nAtoms, params.cloudTemperature, params.cloudRadius,
                                       params.initialPosition, params.initialVelocity);
    // create the detector
    detector = new AISDetector(atomEnsemble, coherenceLength);

    // create the kinematic propagator
    if(params.potentialType == "zero_pot")
    {
        kinematicPropagator = new AISKinematicPropagator(zeroU,zeroGrad,zeroGrad,zeroHess,zeroHess,zeroHess);
    }
    else
    {
        std::cerr << "Potential type "<<params.potentialType<<" not recognized. Exiting." << std::endl;
        exit(1);
    }

    // create the pulse propagators
    using wavefrontFunctionType = double(*)(doubleThreeVector);
    using delWavefrontFunctionType = doubleThreeVector(*)(doubleThreeVector);
    using rabifreqFunctionType  = double(*)(doubleThreeVector, __float128, __float128);

    wavefrontFunctionType wff;
    delWavefrontFunctionType dwff;
    rabifreqFunctionType rff;

    __float128 omega_;
    doubleThreeVector k_;
    double rabiFreq_;

    AISLaserBeam* beam;

    std::cout << "Number of pulses: " << params.rabiFrequencies.size() << std::endl;

    for(int i = 0; i < params.rabiFrequencies.size(); ++i)
    {   
        if(params.wavefrontTypeVector[i] == "flat_square")
        {
            wff = flatWavefront;
            dwff = flatGradientWavefront;
            rff = flatSquareEnvelope;
        }
        else{
            std::cerr << "Wavefront type " << params.wavefrontTypeVector[i] << " not recognized. Exiting." << std::endl;
            exit(1);
        }
        // create the beam
        k_        = {params.kXVector[i], params.kYVector[i], params.kZVector[i]};
        omega_    = params.omegaVector[i];
        rabiFreq_ = params.rabiFrequencies[i];

        beam = new AISLaserBeam(k_,omega_,rabiFreq_,wff,dwff,rff);

        // create the propagator
        AISPulsePropagator* pulsePropagator = new AISPulsePropagator(beam, params.initialPulseTimes[i], params.finalPulseTimes[i], kinematicPropagator);

        // add the propagator to the list
        pulsePropagators.push_back(pulsePropagator);
    }

    // store the params
    this->params = params;
}

AISDriver::~AISDriver()
{
    delete atomEnsemble;
    delete detector;
    delete kinematicPropagator;
    for(auto pulsePropagator : pulsePropagators)
    {
        delete pulsePropagator;
    }
    for(auto laserBeam : laserBeams)
    {
        delete laserBeam;
    }
}

void AISDriver::Run()
{

    for(int i = 0; i < pulsePropagators.size(); ++i)
    {
        // kinematic propagation
        if(i==0 && params.initialPulseTimes[i] != 0.0q) // propagate the ensemble to the initial pulse time
        {
            kinematicPropagator->PropagateEnsemble(atomEnsemble, params.initialPulseTimes[i]);
        }
        // pulse propagation
        pulsePropagators[i]->PropagateEnsemble(atomEnsemble);
    }

    // final propagation
    if(params.finalPropagationTime != params.finalPulseTimes.back())
    {
        kinematicPropagator->PropagateEnsemble(atomEnsemble, params.finalPropagationTime);
    }
}

void AISDriver::Detect()
{
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
            AISPortFrame* portFrame = detector->GetPortFrame(portFrameIndex);
            // get the port
            AISPort* port = portFrame->GetPort(sampledPortIndex);
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
