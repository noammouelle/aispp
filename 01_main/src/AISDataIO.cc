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

std::vector<std::string> readLinesFromFile(std::string fName)
{
    std::ifstream file(fName);
    if (!file.is_open()) {
        std::cerr << "File " << fName << " not found. Returning empty vector" << std::endl;
        return {};
    }

    std::string line, combinedLine;
    std::vector<std::string> lines;

    while (std::getline(file, line)) {
        // Skip comment lines
        if (line.empty() || line[0] == '#') {
            continue;
        }

        // Combine lines that end with a backslash
        if (line.back() == '\\') {
            combinedLine += line.substr(0, line.size() - 1);
        } else {
            combinedLine += line;
            lines.push_back(combinedLine);
            combinedLine.clear();
        }
    }

    file.close();

    return lines;    
}

AISParams readParamsFromFile(std::string fName)
{
    std::ifstream file(fName);
    if (!file.is_open()) {
        std::cerr << "File " << fName << " not found. Returning default params" << std::endl;
        return AISParams();
    }

    // get the lines
    std::vector<std::string> lines = readLinesFromFile(fName);

    // Define the parameter maps
    std::map<std::string, bool> BoolParams;
    std::map<std::string, int> IntParams;
    std::map<std::string, double> DoubleParams;
    std::map<std::string, __float128> QuadParams;
    std::map<std::string, std::string> StrParams;
    std::map<std::string, std::vector<double>> DoubleArrayParams;
    std::map<std::string, std::vector<__float128>> QuadArrayParams;
    std::map<std::string, std::vector<std::string>> StrArrayParams;

    // Define the mandatory parameter keys
    std::set<std::string> MandatoryBoolParamsKeys = {"printprobs","usemcbranching"};
    std::set<std::string> MandatoryIntParamsKeys = {"natoms"};
    std::set<std::string> MandatoryDoubleParamsKeys = {"sigma", "temp", "amplitudethreshold", "coherencelength"};
    std::set<std::string> MandatoryQuadParamsKeys = {"detectiontime"};
    std::set<std::string> MandatoryStrParamsKeys = {"utype"};
    std::set<std::string> MandatoryDoubleArrayParamsKeys = {"x0", "v0", "kx", "ky", "kz","rabifreq", "phi0"};
    std::set<std::string> MandatoryQuadArrayParamsKeys = {"t0", "t1", "omega"};
    std::set<std::string> MandatoryStrArrayParamsKeys = {"wtype"};

    // Define the optional keys
    std::set<std::string> OptBoolParamsKeys = {};
    std::set<std::string> OptIntParamsKeys = {};
    std::set<std::string> OptDoubleParamsKeys = {};
    std::set<std::string> OptQuadParamsKeys = {};
    std::set<std::string> OptStrParamsKeys = {};
    std::set<std::string> OptDoubleArrayParamsKeys = {};
    std::set<std::string> OptQuadArrayParamsKeys = {"beta"};
    std::set<std::string> OptStrArrayParamsKeys = {};

    // Define the possible parameters combining the mandatory and optional keys
    std::set<std::string> BoolParamsKeys;
    std::set<std::string> IntParamsKeys;
    std::set<std::string> DoubleParamsKeys;
    std::set<std::string> QuadParamsKeys;
    std::set<std::string> StrParamsKeys;
    std::set<std::string> DoubleArrayParamsKeys;
    std::set<std::string> QuadArrayParamsKeys;
    std::set<std::string> StrArrayParamsKeys;

    std::set_union(MandatoryBoolParamsKeys.begin(), MandatoryBoolParamsKeys.end(), OptBoolParamsKeys.begin(), OptBoolParamsKeys.end(), std::inserter(BoolParamsKeys, BoolParamsKeys.begin()));
    std::set_union(MandatoryIntParamsKeys.begin(), MandatoryIntParamsKeys.end(), OptIntParamsKeys.begin(), OptIntParamsKeys.end(), std::inserter(IntParamsKeys, IntParamsKeys.begin()));
    std::set_union(MandatoryDoubleParamsKeys.begin(), MandatoryDoubleParamsKeys.end(), OptDoubleParamsKeys.begin(), OptDoubleParamsKeys.end(), std::inserter(DoubleParamsKeys, DoubleParamsKeys.begin()));
    std::set_union(MandatoryQuadParamsKeys.begin(), MandatoryQuadParamsKeys.end(), OptQuadParamsKeys.begin(), OptQuadParamsKeys.end(), std::inserter(QuadParamsKeys, QuadParamsKeys.begin()));
    std::set_union(MandatoryStrParamsKeys.begin(), MandatoryStrParamsKeys.end(), OptStrParamsKeys.begin(), OptStrParamsKeys.end(), std::inserter(StrParamsKeys, StrParamsKeys.begin()));
    std::set_union(MandatoryDoubleArrayParamsKeys.begin(), MandatoryDoubleArrayParamsKeys.end(), OptDoubleArrayParamsKeys.begin(), OptDoubleArrayParamsKeys.end(), std::inserter(DoubleArrayParamsKeys, DoubleArrayParamsKeys.begin()));
    std::set_union(MandatoryQuadArrayParamsKeys.begin(), MandatoryQuadArrayParamsKeys.end(), OptQuadArrayParamsKeys.begin(), OptQuadArrayParamsKeys.end(), std::inserter(QuadArrayParamsKeys, QuadArrayParamsKeys.begin()));
    std::set_union(MandatoryStrArrayParamsKeys.begin(), MandatoryStrArrayParamsKeys.end(), OptStrArrayParamsKeys.begin(), OptStrArrayParamsKeys.end(), std::inserter(StrArrayParamsKeys, StrArrayParamsKeys.begin()));

    std::string line;

    for (const auto& line : lines) {
        // Parse the line
        std::istringstream iss(line);
        std::string key;
        iss >> key;

        // Check and insert into the appropriate map
        if (BoolParamsKeys.find(key) != BoolParamsKeys.end()) {
            bool value;
            iss >> value;
            BoolParams[key] = value;
        }
        else if (IntParamsKeys.find(key) != IntParamsKeys.end()) {
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

    // Check for missing parameters
    for (const auto& key : MandatoryBoolParamsKeys) {
        if (BoolParams.find(key) == BoolParams.end()) {
            throw std::runtime_error("Missing bool parameter: " + key);
        }
    }
    for (const auto& key : MandatoryIntParamsKeys) {
        if (IntParams.find(key) == IntParams.end()) {
            throw std::runtime_error("Missing integer parameter: " + key);
        }
    }
    for (const auto& key : MandatoryDoubleParamsKeys) {
        if (DoubleParams.find(key) == DoubleParams.end()) {
            throw std::runtime_error("Missing double parameter: " + key);
        }
    }
    for (const auto& key : MandatoryQuadParamsKeys) {
        if (QuadParams.find(key) == QuadParams.end()) {
            throw std::runtime_error("Missing quad parameter: " + key);
        }
    }
    for (const auto& key : MandatoryStrParamsKeys) {
        if (StrParams.find(key) == StrParams.end()) {
            throw std::runtime_error("Missing string parameter: " + key);
        }
    }
    for (const auto& key : MandatoryDoubleArrayParamsKeys) {
        if (DoubleArrayParams.find(key) == DoubleArrayParams.end()) {
            throw std::runtime_error("Missing double array parameter: " + key);
        }
    }
    for (const auto& key : MandatoryQuadArrayParamsKeys) {
        if (QuadArrayParams.find(key) == QuadArrayParams.end()) {
            throw std::runtime_error("Missing quad array parameter: " + key);
        }
    }
    for (const auto& key : MandatoryStrArrayParamsKeys) {
        if (StrArrayParams.find(key) == StrArrayParams.end()) {
            throw std::runtime_error("Missing string array parameter: " + key);
        }
    }

    // Sanity checks
    // detection time > finalPulseTimes[-1]
    if (QuadParams["detectiontime"] < QuadArrayParams["t1"].back()) {
        throw std::runtime_error("Detection time must be greater than the final pulse time");
    }

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
    params.betaVector = QuadArrayParams["beta"];

    params.detectionTime = QuadParams["detectiontime"];

    params.amplitudeThreshold = DoubleParams["amplitudethreshold"];
    params.coherenceLength = DoubleParams["coherencelength"];

    params.printPorts = BoolParams["printprobs"];
    params.useMcBranching = BoolParams["usemcbranching"];

    return params;
}