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
    // create the "drift" propagator
    driftPropagator = new AISLinearGravityPropagator(params.interogationTime / params.nSteps);
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
                                                
    // compute the detuned frequency
    __float128 omega1 = computeDetunedOmega(omegaSr87, - 1 * g * params.interogationTime);
    __float128 omega2 = computeDetunedOmega(omegaSr87, - 2 * g * params.interogationTime);
    __float128 omega3 = computeDetunedOmega(omegaSr87, - 3 * g * params.interogationTime);

    // compute the detuned wavevectors
    doubleThreeVector k0 = {0., 0., convertScalarToDouble(omegaSr87 / c)};
    doubleThreeVector k1 = computeDetunedWaveVector(k0, - 1 * g * params.interogationTime);
    doubleThreeVector k2 = computeDetunedWaveVector(k0, - 2 * g * params.interogationTime);
    doubleThreeVector k3 = computeDetunedWaveVector(k0, - 3 * g * params.interogationTime);
    
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

    // set the detector params
    coherenceLength = params.coherenceLength;
}

AISAtomInterferometer::~AISAtomInterferometer(){
    delete atomEnsemble;
    delete beamSplitterPropagator1;
    delete beamSplitterPropagator2;
    delete mirrorPropagator;
    delete driftPropagator;
    delete interrogationTimePropagator;
    delete fpDetector;
}

AISAtomEnsemble* AISAtomInterferometer::GetAtomEnsemble()
{
    return atomEnsemble;
}

void AISAtomInterferometer::run()
{
    // free propagation
    for(int n = 0; n < fParams.nSteps; ++n)
    {
        driftPropagator->PropagateEnsemble(atomEnsemble);
    }
    // beam-splitting
    beamSplitterPropagator1->PropagateEnsemble(atomEnsemble);
    // free propagation
    for(int n = 0; n < fParams.nSteps; ++n)
    {
        interrogationTimePropagator->PropagateEnsemble(atomEnsemble);
    }
    // mirror
    mirrorPropagator->PropagateEnsemble(atomEnsemble);
    // free propagation
    for(int n = 0; n < fParams.nSteps; ++n)
    {
        interrogationTimePropagator->PropagateEnsemble(atomEnsemble);
    }
    // beam-splitting
    beamSplitterPropagator2->PropagateEnsemble(atomEnsemble);
    // free propagation
    for(int n = 0; n < fParams.nSteps; ++n)
    {
    driftPropagator->PropagateEnsemble(atomEnsemble);
    }
}

void AISAtomInterferometer::detect()
{
    // create the detector
    fpDetector = new AISDetector(atomEnsemble, coherenceLength);
    // detect the atoms
    sampledPortIndices = fpDetector->SampleAllPorts();
}

void AISAtomInterferometer::write(std::string filename)
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