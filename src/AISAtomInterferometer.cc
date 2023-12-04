#include "AISAtomInterferometerParams.hh"
#include "AISAtomInterferometer.hh"
#include "AISAtomEnsemble.hh"
#include "AISLinearGravityPropagator.hh"
#include "AISIntensityProfile.hh"
#include "AIStttPropagator.hh"
#include "AISWaveFront.hh"
#include "AISConstants.hh"

AISAtomInterferometer::AISAtomInterferometer(AISAtomInterferometerParams params)
{
    /*
    Initializes a 3 pulse atom interferometer. This will be upgraded for general 
    LMT atom interferometer of order N.
    */

    fParams = params;

    // create the ensemble
    atomEnsemble = new AISAtomEnsemble(params.nAtoms, params.cloudTemperature, params.cloudWidth,
                                       params.initialPosition, params.initialVelocity);
    // create the "drift" propagators
    initialDriftPropagator = new AISLinearGravityPropagator(params.initialPropagationTime / params.nSteps);
    finalDriftPropagator   = new AISLinearGravityPropagator(params.finalPropagationTime / params.nSteps);
    // create the interrogation time propagators
    interrogationTimePropagator = new AISLinearGravityPropagator(params.interogationTime / params.nSteps);
    interrogationTimePropagator->SetAddEnergyPhase(true);
    
    // create the intensity profile
    AISIntensityProfile* intensityProfile = new AISIntensityProfile(params.beamRadius, params.rabiFrequency);
    // create the wavefront
    AISWaveFront* wavefront    = new AISWaveFront(zeroAberrationFunction, zeroAberrationFunction,
                                                  {0.0, 0.0, 0.0}, 0.0);
    AISWaveFront* wavefrontPsr = new AISWaveFront(zeroAberrationFunction, zeroAberrationFunction,
                                                  params.psrGradient, params.laserPhase);

    // compute the vertical velocities at each pulse
    double vzBeamSplitter1 = params.initialVelocity[2] - 1 * g * params.initialPropagationTime; 

    double vzLmt1 = vzBeamSplitter1 - g * params.lmtDelayTime; 
    double vzLmt2 = vzLmt1 - g * ((params.lmtOrder - 1) / 2 * params.lmtDelayTime + params.interogationTime);

    double vzMirror = vzLmt2 - g * (params.lmtOrder - 1) / 2 * params.lmtDelayTime; 

    double vzLmt3 = vzMirror - g * params.lmtDelayTime; 
    double vzLmt4 = vzLmt3 - g * ((params.lmtOrder - 1) / 2 * params.lmtDelayTime + params.interogationTime); 

    double vzBeamSplitter2 = vzLmt4 - g * (params.lmtOrder - 1) / 2 * params.lmtDelayTime;; 
                                                
    // compute the detuned frequency
    __float128 omega1 = computeDetunedOmega(omegaSr87, vzBeamSplitter1);
    __float128 omega2 = computeDetunedOmega(omegaSr87, vzMirror);
    __float128 omega3 = computeDetunedOmega(omegaSr87, vzBeamSplitter2);

    // compute the detuned wavevectors
    doubleThreeVector k0 = {0., 0., convertScalarToDouble(omegaSr87 / c)};
    doubleThreeVector k1 = computeDetunedWaveVector(k0, vzBeamSplitter1);
    doubleThreeVector k2 = computeDetunedWaveVector(k0, vzMirror);
    doubleThreeVector k3 = computeDetunedWaveVector(k0, vzBeamSplitter2);
    
    // create the first pi/2 pulse ttt propagator
    beamSplitterPropagator1 = new AIStttPropagator(pi128 / (2 * params.rabiFrequency));
    beamSplitterPropagator1->SetWaveFront(wavefront);
    beamSplitterPropagator1->SetIntensityProfile(intensityProfile);
    beamSplitterPropagator1->SetOmega(omega1);
    beamSplitterPropagator1->SetWaveVector(k1);
    // create the mirror pulse
    mirrorPropagator = new AIStttPropagator(pi128 / (params.rabiFrequency));
    mirrorPropagator->SetWaveFront(wavefront);
    mirrorPropagator->SetIntensityProfile(intensityProfile);
    mirrorPropagator->SetOmega(omega2);
    mirrorPropagator->SetWaveVector(k2);
    // create the second pi/2 pulse ttt propagator
    beamSplitterPropagator2 = new AIStttPropagator(pi128 / (2 * params.rabiFrequency));
    beamSplitterPropagator2->SetWaveFront(wavefrontPsr);
    beamSplitterPropagator2->SetIntensityProfile(intensityProfile);
    beamSplitterPropagator2->SetOmega(omega3);
    beamSplitterPropagator2->SetWaveVector(k3);  

    // set the LMT propagator params
    lmtPropagator1 = new AISLmtTttPropagator(params.lmtOrder, vzLmt1, params.lmtDelayTime, params.interogationTime, 0);
    lmtPropagator2 = new AISLmtTttPropagator(params.lmtOrder, vzLmt2, params.lmtDelayTime, params.interogationTime, 1);
    lmtPropagator3 = new AISLmtTttPropagator(params.lmtOrder, vzLmt3, params.lmtDelayTime, params.interogationTime, 2);
    lmtPropagator4 = new AISLmtTttPropagator(params.lmtOrder, vzLmt4, params.lmtDelayTime, params.interogationTime, 3);

    lmtPropagator1->SetWaveFronts(wavefront, wavefront);
    lmtPropagator2->SetWaveFronts(wavefront, wavefront);
    lmtPropagator3->SetWaveFronts(wavefront, wavefront);
    lmtPropagator4->SetWaveFronts(wavefront, wavefront);

    lmtPropagator1->SetIntensityProfiles(intensityProfile, intensityProfile);
    lmtPropagator2->SetIntensityProfiles(intensityProfile, intensityProfile);
    lmtPropagator3->SetIntensityProfiles(intensityProfile, intensityProfile);
    lmtPropagator4->SetIntensityProfiles(intensityProfile, intensityProfile);

    // set the detector params
    coherenceLength = params.coherenceLength;
}

AISAtomInterferometer::~AISAtomInterferometer(){
    delete atomEnsemble;
    delete beamSplitterPropagator1;
    delete beamSplitterPropagator2;
    delete mirrorPropagator;
    delete initialDriftPropagator;
    delete finalDriftPropagator;
    delete interrogationTimePropagator;
    delete fpDetector;
    delete lmtPropagator1;
    delete lmtPropagator2;
    delete lmtPropagator3;
    delete lmtPropagator4;
}

AISAtomEnsemble* AISAtomInterferometer::GetAtomEnsemble()
{
    return atomEnsemble;
}

void AISAtomInterferometer::run()
{
    // loop over the atoms
    #pragma omp parallel for
    for(int atomIndex = 0; atomIndex < atomEnsemble->GetNumberOfAtoms(); ++atomIndex)
    {
        AISAtom* currentAtom = atomEnsemble->GetAtom(atomIndex);
        runAtom(currentAtom);
    }
}

void AISAtomInterferometer::runAtom(AISAtom* anAtom)
{
    // free propagation
    for(int n = 0; n < fParams.nSteps; ++n)
    {
        initialDriftPropagator->PropagateAtom(anAtom);
    }

    // beam-splitting
    beamSplitterPropagator1->PropagateAtom(anAtom);

    // first LMT block
    lmtPropagator1->PropagateAtom(anAtom);

    // free propagation
    for(int n = 0; n < fParams.nSteps; ++n)
    {
        interrogationTimePropagator->PropagateAtom(anAtom);
    }

    // second LMT block
    lmtPropagator2->PropagateAtom(anAtom);

    // mirror
    mirrorPropagator->PropagateAtom(anAtom);

    // third LMT block
    lmtPropagator3->PropagateAtom(anAtom);

    // free propagation
    for(int n = 0; n < fParams.nSteps; ++n)
    {
        interrogationTimePropagator->PropagateAtom(anAtom);
    }

    // fourth LMT block
    lmtPropagator4->PropagateAtom(anAtom);

    // beam-splitting
    beamSplitterPropagator2->PropagateAtom(anAtom);
    // free propagation
    for(int n = 0; n < fParams.nSteps; ++n)
    {
    finalDriftPropagator->PropagateAtom(anAtom);
    }
}

void AISAtomInterferometer::detect()
{
    // create the detector
    fpDetector = new AISDetector(atomEnsemble, coherenceLength);
    // detect the atoms
    sampledPortIndices = fpDetector->SampleAllPorts();
}

void AISAtomInterferometer::writeDetectedAtomsInfo(std::string filename)
{
    // open a h5 file
    H5::H5File file(filename, H5F_ACC_TRUNC);
    int numSamples = fpDetector->GetNumberOfSamples();

    // create the data vectors
    std::vector<doubleThreeVector> positions(numSamples);
    std::vector<doubleThreeVector> velocities(numSamples);
    std::vector<int> states(numSamples);
    std::vector<double> phaseShifts(numSamples);

    // loop over the samples to fill in the data
    int sampleIndex = 0;
    for(int portFrameIndex = 0; portFrameIndex < fpDetector->GetNumberOfPortFrames(); portFrameIndex++)
    {    
        int sampledPortIndex = sampledPortIndices[portFrameIndex];
        if(sampledPortIndex != -1)
        {
            // get the port frame
            AISPortFrame* portFrame = fpDetector->GetPortFrame(portFrameIndex);
            // get the port
            AISPort* port = portFrame->GetPort(sampledPortIndex);
            // get the data
            doubleThreeVector currentPosition = port->position;
            doubleThreeVector currentVelocity = port->velocity;
            int currentState = port->state;
            double currentPhaseShift = port->phaseShift;
            // fill in the data
            for(int i = 0; i < 3; ++i)
            {
                positions[sampleIndex][i] = currentPosition[i];
                velocities[sampleIndex][i] = currentVelocity[i];
            }
            states[sampleIndex] = currentState;
            phaseShifts[sampleIndex] = currentPhaseShift;

            sampleIndex++;
        }        
    }
    // define the dimensionality of the datasets
    hsize_t numSamplesHsize = numSamples;
    hsize_t dim_positions[2] = {numSamplesHsize, 3};
    hsize_t dim_velocities[2] = {numSamplesHsize, 3};
    hsize_t dim_states[1] = {numSamplesHsize};
    hsize_t dim_phaseShifts[1] = {numSamplesHsize};
    // create the dataspace
    H5::DataSpace dataspace_positions(2, dim_positions);
    H5::DataSpace dataspace_velocities(2, dim_velocities);
    H5::DataSpace dataspace_states(1, dim_states);
    H5::DataSpace dataspace_phaseShifts(1, dim_phaseShifts);
    // create the datasets
    H5::DataSet dataset_positions = file.createDataSet("positions", H5::PredType::NATIVE_DOUBLE, dataspace_positions);
    H5::DataSet dataset_velocities = file.createDataSet("velocities", H5::PredType::NATIVE_DOUBLE, dataspace_velocities);
    H5::DataSet dataset_states = file.createDataSet("states", H5::PredType::NATIVE_INT, dataspace_states);
    H5::DataSet dataset_phaseShifts = file.createDataSet("phaseShifts", H5::PredType::NATIVE_DOUBLE, dataspace_phaseShifts);
    // write the data
    dataset_positions.write(positions.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_velocities.write(velocities.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_states.write(states.data(), H5::PredType::NATIVE_INT);
    dataset_phaseShifts.write(phaseShifts.data(), H5::PredType::NATIVE_DOUBLE);

    // close the file
    file.close();
}

void AISAtomInterferometer::writeWavepacketsInfo(std::string filename)
{
    /* This was written only for the case where each atom has the 
     * same number of wavepackets. (i.e. no wavepackets are deleted)
    */

    // open H5 file
    H5::H5File file(filename, H5F_ACC_TRUNC);

    // get the number of atoms
    int nAtoms = atomEnsemble->GetNumberOfAtoms();
    // get the number of wavepackets
    int nWavepackets = atomEnsemble->GetAtom(0)->GetNumberOfWavePackets();

    // create data vectors
    std::vector<std::vector<doubleThreeVector>> positions(nAtoms);
    std::vector<std::vector<doubleThreeVector>> velocities(nAtoms);
    std::vector<std::vector<double>> phases(nAtoms);
    std::vector<std::vector<__float128>> quadPhases(nAtoms);
    std::vector<std::vector<double>> amplitudes(nAtoms);
    std::vector<std::vector<int>> states(nAtoms);

    // loop over the atoms
    for(int atomIndex = 0; atomIndex < nAtoms; atomIndex++){
        // get the atom
        AISAtom* currentAtom = atomEnsemble->GetAtom(atomIndex);

        // loop over the wavepacket vector
        for(int wavepacketIndex = 0; wavepacketIndex < nWavepackets; wavepacketIndex++){
            // get the wavepacket
            AISWavePacket* currentWavepacket = currentAtom->GetWavePacket(wavepacketIndex);

            // get the data
            doubleThreeVector currentPosition = currentWavepacket->GetPosition();
            doubleThreeVector currentVelocity = currentWavepacket->GetVelocity();
            double currentPhase = currentWavepacket->GetPhaseDouble();
            __float128 currentQuadPhase = currentWavepacket->GetPhaseQuad();
            double currentAmplitude = currentWavepacket->GetAmplitude();
            int currentState = currentWavepacket->GetState();

            // fill in the data
            positions[atomIndex].push_back(currentPosition);
            velocities[atomIndex].push_back(currentVelocity);
            phases[atomIndex].push_back(currentPhase);
            quadPhases[atomIndex].push_back(currentQuadPhase);
            amplitudes[atomIndex].push_back(currentAmplitude);
            states[atomIndex].push_back(currentState);
        }
    }

    // define the dimensionality of the datasets
    hsize_t nAtomsHsize       = nAtoms;
    hsize_t nWavepacketsHsize = nWavepackets;
    hsize_t dim_positions[3]  = {nAtomsHsize, nWavepacketsHsize, 3};
    hsize_t dim_velocities[3] = {nAtomsHsize, nWavepacketsHsize, 3};
    hsize_t dim_phases[2]     = {nAtomsHsize, nWavepacketsHsize};
    hsize_t dim_quadPhases[2] = {nAtomsHsize, nWavepacketsHsize};
    hsize_t dim_amplitudes[2] = {nAtomsHsize, nWavepacketsHsize};
    hsize_t dim_states[2]     = {nAtomsHsize, nWavepacketsHsize};
    
    // create the dataspace
    H5::DataSpace dataspace_positions(3, dim_positions);
    H5::DataSpace dataspace_velocities(3, dim_velocities);
    H5::DataSpace dataspace_phases(2, dim_phases);
    H5::DataSpace dataspace_quadPhases(2, dim_quadPhases);
    H5::DataSpace dataspace_amplitudes(2, dim_amplitudes);
    H5::DataSpace dataspace_states(2, dim_states);
    
    // create the datasets
    H5::DataSet dataset_positions = file.createDataSet("positions", H5::PredType::NATIVE_DOUBLE, dataspace_positions);
    H5::DataSet dataset_velocities = file.createDataSet("velocities", H5::PredType::NATIVE_DOUBLE, dataspace_velocities);
    H5::DataSet dataset_phases = file.createDataSet("phases", H5::PredType::NATIVE_DOUBLE, dataspace_phases);
    H5::DataSet dataset_quadPhases = file.createDataSet("quadPhases", H5::PredType::NATIVE_LDOUBLE, dataspace_quadPhases);
    H5::DataSet dataset_amplitudes = file.createDataSet("amplitudes", H5::PredType::NATIVE_DOUBLE, dataspace_amplitudes);
    H5::DataSet dataset_states = file.createDataSet("states", H5::PredType::NATIVE_INT, dataspace_states);

    // write the data
    dataset_positions.write(positions.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_velocities.write(velocities.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_phases.write(phases.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_quadPhases.write(quadPhases.data(), H5::PredType::NATIVE_LDOUBLE);
    dataset_amplitudes.write(amplitudes.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_states.write(states.data(), H5::PredType::NATIVE_INT);

    // close the file
    file.close();
}

doubleThreeVector AISAtomInterferometer::computeDetunedWaveVector(const doubleThreeVector& k0, const double& vz)
{
    doubleThreeVector detunedK = k0;
    for(int i = 0; i < 3; ++i)
    {
        if(i == 2){
            detunedK[i] = k0[i] /(1 - vz / c);
        }
    }

    return detunedK;
}

__float128 AISAtomInterferometer::computeDetunedOmega(const __float128& omega0, const double& vz){
    return omega0 / (1 - vz / c);
}