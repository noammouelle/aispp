#include "AISDataIO.hh"

void writeAtomEnsembleToFile(std::string fName, std::unique_ptr<AISAtomEnsemble>& atomEnsemble)
{
    // write a file to cross check
    std::ofstream outFile(fName);
    std::cout << "Writing to file: " << fName << std::endl;
    outFile << std::fixed; // Use fixed-point notation
    outFile << std::setprecision(15);

    outFile << "State, Amplitude, X, Y, Z, VX, VY, VZ, Phase, PhaseErr, PhaseQuad, Time" << std::endl;
    for(int i = 0; i < atomEnsemble->GetNumberOfAtoms(); ++i)
    {
        std::unique_ptr<AISAtom>& currentAtom = atomEnsemble->GetAtom(i);
        for(int j = 0; j < currentAtom->GetNumberOfWavePackets(); ++j)
        {
            std::unique_ptr<AISWavePacket>& currentWavePacket = currentAtom->GetWavePacket(j);
            int currentState = currentWavePacket->GetState();
            double currentAmplitude = currentWavePacket->GetAmplitude();
            doubleThreeVector currentPos = currentWavePacket->GetPosition();
            doubleThreeVector currentVel = currentWavePacket->GetVelocity();
            double currentPhaseDouble = currentWavePacket->GetPhaseDouble();
            double currentPhaseDoubleErr = currentWavePacket->GetPhaseDoubleError();
            __float128 currentPhaseQuad = currentWavePacket->GetPhaseQuad();
            __float128 currentTime = currentWavePacket->GetTime();

            long double currentPhaseQuadAsLongDouble = static_cast<long double>(currentPhaseQuad);
            long double currentTimeAsLongDouble = static_cast<long double>(currentTime);

            char currentPhaseQuadStr[50];  // Adjust the buffer size as needed
            char currentTimeStr[50];  // Adjust the buffer size as needed
            snprintf(currentPhaseQuadStr, sizeof(currentPhaseQuadStr), "%.34Le", currentPhaseQuadAsLongDouble);
            snprintf(currentTimeStr, sizeof(currentTimeStr), "%.34Le", currentTimeAsLongDouble);

            outFile << currentState << "," << currentAmplitude << "," << currentPos[0] << "," << currentPos[1] << "," << currentPos[2] << "," 
                    << currentVel[0] << "," << currentVel[1] << "," << currentVel[2] << ","
                    << currentPhaseDouble << "," << currentPhaseDoubleErr << "," << currentPhaseQuadStr << "," << currentTimeStr << std::endl;
        }
    }

    outFile.close();
}

AISParams readParamsFromFile(std::string fName)
{
    std::ifstream file(fName);
    if (!file.is_open()) {
        std::cerr << "File " << fName << " not found. Returning default params" << std::endl;
        return AISParams();
    }

    // Define the parameter maps
    std::map<std::string, int> IntParams;
    std::map<std::string, double> DoubleParams;
    std::map<std::string, __float128> QuadParams;
    std::map<std::string, std::string> StrParams;
    std::map<std::string, std::vector<double>> DoubleArrayParams;
    std::map<std::string, std::vector<__float128>> QuadArrayParams;
    std::map<std::string, std::vector<std::string>> StrArrayParams;

    // Define the possible parameter keys using sets
    std::set<std::string> IntParamsKeys = {"natoms"};
    std::set<std::string> DoubleParamsKeys = {"sigma", "temp"};
    std::set<std::string> QuadParamsKeys = {"finalproptime"};
    std::set<std::string> StrParamsKeys = {"utype"};
    std::set<std::string> DoubleArrayParamsKeys = {"x0", "v0", "kx", "ky", "kz","rabifreq", "phi0"};
    std::set<std::string> QuadArrayParamsKeys = {"t0", "t1", "omega"};
    std::set<std::string> StrArrayParamsKeys = {"wtype"};

    std::string line;

    while (std::getline(file, line)) {
        // Skip comment lines
        if (line.empty() || line[0] == '#') {
            continue;
        }

        std::istringstream iss(line);
        std::string key;
        iss >> key;

        // Check and insert into the appropriate map
        if (IntParamsKeys.find(key) != IntParamsKeys.end()) {
            int value;
            iss >> value;
            IntParams[key] = value;
        }
        else if (DoubleParamsKeys.find(key) != DoubleParamsKeys.end()) {
            double value;
            iss >> value;
            DoubleParams[key] = value;
        } else if (QuadParamsKeys.find(key) != QuadParamsKeys.end()) {
            __float128 value;
            std::string valueStr;
            iss >> valueStr;
            value = strtoflt128(valueStr.c_str(), NULL);
            QuadParams[key] = value;
        } else if (StrParamsKeys.find(key) != StrParamsKeys.end()) {
            std::string value;
            iss >> value;
            StrParams[key] = value;
        } else if (DoubleArrayParamsKeys.find(key) != DoubleArrayParamsKeys.end()) {
            std::vector<double> values;
            double value;
            while (iss >> value) {
                values.push_back(value);
            }
            DoubleArrayParams[key] = values;
        } else if (QuadArrayParamsKeys.find(key) != QuadArrayParamsKeys.end()) {
            std::vector<__float128> values;
            std::string valueStr;
            while (iss >> valueStr) {
                values.push_back(strtoflt128(valueStr.c_str(), NULL));
            }
            QuadArrayParams[key] = values;
        } else if (StrArrayParamsKeys.find(key) != StrArrayParamsKeys.end()) {
            std::vector<std::string> values;
            std::string value;
            while (iss >> value) {
                values.push_back(value);
            }
            StrArrayParams[key] = values;
        }
        else {
            std::cerr << "Unknown key: " << key << " Aborting." << std::endl;
            exit(1);
        }
    }

    file.close();

    AISParams params;
    params.nAtoms           = IntParams["natoms"];
    params.cloudRadius      = DoubleParams["sigma"];
    params.cloudTemperature = DoubleParams["temp"];
    params.potentialType    = StrParams["utype"];
    params.initialPosition  = {DoubleArrayParams["x0"][0], DoubleArrayParams["x0"][1], DoubleArrayParams["x0"][2]};
    params.initialVelocity  = {DoubleArrayParams["v0"][0], DoubleArrayParams["v0"][1], DoubleArrayParams["v0"][2]};
    params.initialPulseTimes = QuadArrayParams["t0"];
    params.finalPulseTimes = QuadArrayParams["t1"];
    params.kXVector = DoubleArrayParams["kx"];
    params.kYVector = DoubleArrayParams["ky"];
    params.kZVector = DoubleArrayParams["kz"];

    // Rabi frequencies (convert Hz to rad/s)
    for(int i=0; i<DoubleArrayParams["rabifreq"].size(); i++)
    {
        DoubleArrayParams["rabifreq"][i] = DoubleArrayParams["rabifreq"][i] * 2 * pi;
    }
    params.rabiFrequencies = DoubleArrayParams["rabifreq"];
    
    params.phi0 = DoubleArrayParams["phi0"];
    params.omegaVector = QuadArrayParams["omega"];
    params.wavefrontTypeVector = StrArrayParams["wtype"];

    params.finalPropagationTime = QuadParams["finalproptime"];

    return params;
}