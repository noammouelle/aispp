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



}