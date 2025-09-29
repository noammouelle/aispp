#include "AISLaserBeam.hh"

AISLaserBeam::AISLaserBeam(doubleThreeVector aK, __float128 aOmega, double aRabiFreq, double phi0) : k(aK), omega(aOmega), rabiFreq(aRabiFreq), phi0(phi0)
{}

AISLaserBeam::~AISLaserBeam()
{}

doubleThreeVector AISLaserBeam::GetK(__float128 t)
{
    return matrixAdd(this->k,scalarMultiply(this->kChirp,0.5 * static_cast<double>(t)));
}

void AISLaserBeam::SetK(doubleThreeVector kValue)
{
    k = kValue;
}

doubleThreeVector AISLaserBeam::GetKChirp()
{
    return kChirp;
}
void AISLaserBeam::SetKChirp(doubleThreeVector kChirpValue)
{
    kChirp = kChirpValue;
}

/**
 * @brief Calculates the laser frequency at a given time.
 *
 * This function computes the angular frequency (omega) at a specific time (t)
 * by adding the initial angular frequency (omega) to the product of the 
 * frequency chirp rate and the time.
 *
 * @param t The time at which to calculate the angular frequency.
 * @return The angular frequency at time t.
 */
__float128 AISLaserBeam::GetOmega(__float128 t)
{
    return omega + 0.5q * frequencyChirp * t;
}

void AISLaserBeam::SetOmega(__float128 omegaValue)
{
    this->omega=omegaValue;
}

__float128 AISLaserBeam::GetFrequencyChirp()
{
    return frequencyChirp;
}
void AISLaserBeam::SetFrequencyChirp(__float128 frequencyChirpValue)
{
    frequencyChirp=frequencyChirpValue;
}

double AISLaserBeam::GetPhi(const doubleThreeVector& pos)
{
    // check the beam type
    if(beamType == "flat_square")
    {
        return phi0;
    }
    else if(beamType == "gaussian")
    {       
        double spatiallyVaryingPhase = trilinearInterpolation(pos[0], pos[1], pos[2], 
            interpolationGridX, interpolationGridY, interpolationGridZ,
            phaseInterpolationGridValues, interpolationGridX.size(), interpolationGridY.size(), interpolationGridZ.size());

        return phi0 + spatiallyVaryingPhase;
    }
}

double AISLaserBeam::GetPhiWrapper(double xi, void* params)
{
    getPhiWrapperParams* p = (getPhiWrapperParams*) params;
    
    doubleThreeVector pos = p->pos;
    int index = p->index;
    std::shared_ptr<AISLaserBeam> beam = p->beam;
    // only modify the selected index
    pos[index] = xi;
    return beam->GetPhi(pos);
}

doubleThreeVector AISLaserBeam::GetDelPhi(const doubleThreeVector& pos)
{
    return {0., 0., 0.}; // TODO: implement with numerical differentiation
}

double AISLaserBeam::GetRabiFreq(const doubleThreeVector& pos, const __float128& t0, const __float128& t)
{
    // check the beam type
    double effectiveRabiFreq = 0.0;
    if(beamType == "flat_square")
    {
        effectiveRabiFreq = rabiFreq; // constant rabi freq
    }
    else if(beamType == "gaussian")
    {
        double amplitude = trilinearInterpolation(pos[0], pos[1], pos[2], 
            interpolationGridX, interpolationGridY, interpolationGridZ, 
            amplitudeInterpolationGridValues, interpolationGridX.size(), interpolationGridY.size(), interpolationGridZ.size());
        
        return amplitude * rabiFreq;
    }
    return effectiveRabiFreq;
}

double AISLaserBeam::GetW0()
{
    return w0;
}
void AISLaserBeam::SetW0(double w0Value)
{
    w0 = w0Value;
}

double AISLaserBeam::GetZLaser()
{
    return zLaser;
}
void AISLaserBeam::SetZLaser(double zLaserValue)
{
    zLaser = zLaserValue;
}

double AISLaserBeam::GetFocalLength()
{
    return focalLength;
}
void AISLaserBeam::SetFocalLength(double focalLengthValue)
{
    focalLength = focalLengthValue;
}

void AISLaserBeam::SetZernikeCoeffs(std::map<int, double> zernikeCoeffsValue)
{
    zernikeCoeffs = zernikeCoeffsValue;
}
std::map<int, double> AISLaserBeam::GetZernikeCoeffs()
{
    return zernikeCoeffs;
}

void AISLaserBeam::SetBeamRadius(double beamRadiusValue)
{
    beamRadius = beamRadiusValue;
}
double AISLaserBeam::GetBeamRadius()
{
    return beamRadius;
}

void AISLaserBeam::SetBeamType(std::string beamTypeValue)
{
    // check if 'flat_square' or 'gaussian'
    if(beamTypeValue == "flat_square" || beamTypeValue == "gaussian")
    {
        beamType = beamTypeValue;
    }
    else
    {
        std::cerr << "Beam type " << beamTypeValue << " not recognized. Exiting." << std::endl;
        exit(1);
    }
}
std::string AISLaserBeam::GetBeamType()
{
    return beamType;
}

void AISLaserBeam::SetBaselineLength(double baselineLengthValue)
{
    baselineLength = baselineLengthValue;
}
double AISLaserBeam::GetBaselineLength()
{
    return baselineLength;
}

void AISLaserBeam::SetInterpolationGrids(std::string filename)
{
    // Open the HDF5 file
    H5::H5File file(filename, H5F_ACC_RDONLY);

    // Read the interpolation grids
    H5::DataSet datasetX = file.openDataSet("x");
    H5::DataSet datasetY = file.openDataSet("y");
    H5::DataSet datasetZ = file.openDataSet("z");
    H5::DataSet datasetPhase = file.openDataSet("phase");
    H5::DataSet datasetAmplitude = file.openDataSet("amplitude");

    // Get the dimensions of each dataset
    H5::DataSpace dataspaceX = datasetX.getSpace();
    H5::DataSpace dataspaceY = datasetY.getSpace();
    H5::DataSpace dataspaceZ = datasetZ.getSpace();
    H5::DataSpace dataspacePhase = datasetPhase.getSpace();
    H5::DataSpace dataspaceAmplitude = datasetAmplitude.getSpace();

    hsize_t dimsX[1], dimsY[1], dimsZ[1], dimsPhase[1], dimsAmplitude[1];
    dataspaceX.getSimpleExtentDims(dimsX);
    dataspaceY.getSimpleExtentDims(dimsY);
    dataspaceZ.getSimpleExtentDims(dimsZ);
    dataspacePhase.getSimpleExtentDims(dimsPhase);
    dataspaceAmplitude.getSimpleExtentDims(dimsAmplitude);

    // Resize the vectors to match the dataset dimensions
    interpolationGridX.resize(dimsX[0]);
    interpolationGridY.resize(dimsY[0]);
    interpolationGridZ.resize(dimsZ[0]);
    phaseInterpolationGridValues.resize(dimsPhase[0]);
    amplitudeInterpolationGridValues.resize(dimsAmplitude[0]);

    // Read the data into the vectors
    datasetX.read(interpolationGridX.data(), H5::PredType::NATIVE_DOUBLE);
    datasetY.read(interpolationGridY.data(), H5::PredType::NATIVE_DOUBLE);
    datasetZ.read(interpolationGridZ.data(), H5::PredType::NATIVE_DOUBLE);
    datasetPhase.read(phaseInterpolationGridValues.data(), H5::PredType::NATIVE_DOUBLE);
    datasetAmplitude.read(amplitudeInterpolationGridValues.data(), H5::PredType::NATIVE_DOUBLE);
    // Close the file
    file.close();
}