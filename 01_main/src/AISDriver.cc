#include "AISDriver.hh"

AISDriver::AISDriver(AISParams params)
{   
    // three times the De Broglie wavelength of the cloud
    double coherenceLength = 3 * sqrt(2 * pi) * hbar / sqrt(massSr87 * kB * params.cloudTemperature);
    // create the atom ensemble
    atomEnsemble = new AISAtomEnsemble(params.nAtoms, params.cloudRadius, params.cloudTemperature,
                                       params.initialPosition, params.initialVelocity);
    // create the detector
    detector = new AISDetector(atomEnsemble, coherenceLength);

    // create the kinematic propagator
    if(this->params.potentialType == "zero_pot")
    {
        kinematicPropagator = new AISKinematicPropagator(zeroU,zeroGrad,zeroGrad,zeroHess,zeroHess,zeroHess);
    }
    else
    {
        std::cerr << 'Potential type "'<<params.potentialType<<'" not recognized. Exiting.' << std::endl;
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
        AISPulsePropagator* pulsePropagator = new AISPulsePropagator(beam, params.initialTimes[i], params.finalTimes[i], kinematicPropagator);

        // add the propagator to the list
        pulsePropagators.push_back(pulsePropagator);
    }
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

AISDriver::Run()
{
    for(int i = 0; i < pulsePropagators.size(); ++i)
    {
        kinematicPropagator->PropagateEnsemble(atomEnsemble, params.finalTimes[i]);
        pulsePropagators[i]->PropagateEnsemble(atomEnsemble);        
    }
}

AISDriver::Detect()
{
    detector->SampleAllPorts();
}

AISDriver::WriteDetectedAtomsToFile(std::string fName)
{
    
}

AISDriver::WriteWavePacketsToFile(std::string fName)
{
    writeAtomEnsembleToFile(fName, atomEnsemble);
}
