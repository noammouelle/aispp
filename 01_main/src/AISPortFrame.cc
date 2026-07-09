#include "AISPortFrame.hh"

AISPortFrame::AISPortFrame(std::unique_ptr<AISAtom>& anAtom, double coherenceLength, bool ultraFast)
    : ultraFast(ultraFast)
{
    // initialize the port vector
    fpPortVector = std::make_unique<portVector>();

    // for each wavepacket, find the indices of the wavepackets that are close enough to interfere
    std::vector<std::vector<int>> wavePacketIndexGroups;
    for(int wavepacketIndex = 0; wavepacketIndex < anAtom->GetNumberOfWavePackets(); wavepacketIndex++)
    {   
        std::vector<int> wavepacketIndexGroup = createGroup(anAtom, wavepacketIndex, coherenceLength);
        wavePacketIndexGroups.push_back(wavepacketIndexGroup);
    }

    // if two groups are the same, delete one of them
    std::vector<std::vector<int>> uniqueWavePacketIndexGroups;
    for(int i = 0; i < wavePacketIndexGroups.size(); i++)
    {
        if(std::find(uniqueWavePacketIndexGroups.begin(), uniqueWavePacketIndexGroups.end(), wavePacketIndexGroups[i]) == uniqueWavePacketIndexGroups.end())
        {
            uniqueWavePacketIndexGroups.push_back(wavePacketIndexGroups[i]);
        }
    }

    // get the indices of the wavepackets that are close enough to interfere
    std::vector<int> adjacentWavepacketIndices;
    for(int i = 0; i < uniqueWavePacketIndexGroups.size(); i++)
    {
        if (uniqueWavePacketIndexGroups[i].size() > 1)
        {
            for(int j = 0; j < uniqueWavePacketIndexGroups[i].size(); j++)
            {
                adjacentWavepacketIndices.push_back(uniqueWavePacketIndexGroups[i][j]);
            }
        }
    }
    
    // get the paths of the interfering atoms
    for(int i = 0; i < adjacentWavepacketIndices.size(); i++)
    {
        std::unique_ptr<AISWavePacket>& wavePacket = anAtom->GetWavePacket(adjacentWavepacketIndices[i]);
        interferingPaths.push_back(wavePacket->GetPath());
    }

    // create the port vector
    int numberOfPorts = uniqueWavePacketIndexGroups.size();
    //if(numberOfPorts < 1){
        //std::cout << "No ports found for current atom " << std::endl;
        //exit(1);
    //}
    for(int groupIndex = 0; groupIndex < uniqueWavePacketIndexGroups.size(); groupIndex++)
    {
        // create a new port
        std::unique_ptr<AISPort> newPort(new AISPort());
        for(int wavepacketIndex : uniqueWavePacketIndexGroups[groupIndex])
        {
            newPort->wavePacketIndices.push_back(wavepacketIndex);
        }
        fpPortVector->push_back(std::move(newPort));
    }

    // set the port parameters
    setPortParameters(anAtom);
}

AISPortFrame::~AISPortFrame()
{
    // delete all the unique pointers to ports then free the vector smart pointer
    fpPortVector->clear();
    fpPortVector.reset();
}

std::vector<int> AISPortFrame::createGroup(std::unique_ptr<AISAtom>& anAtom, int wavePacketIndex, double coherenceLength)
{
    // create a list of wavepacket indices for current group
    std::vector<int> wavepacketIndexGroup;

    // get the current wavepacket
    wavepacketIndexGroup.push_back(wavePacketIndex);

    // get the position of the current wavepacket
    std::unique_ptr<AISWavePacket>& currentWavePacket = anAtom->GetWavePacket(wavePacketIndex);
    doubleThreeVector currentPosition = currentWavePacket->GetPosition();

    // loop over all the other wavepackets
    for(int otherWavepacketIndex = 0; otherWavepacketIndex < anAtom->GetNumberOfWavePackets(); otherWavepacketIndex++)
    {
        // get the other wavepacket
        std::unique_ptr<AISWavePacket>& otherWavePacket = anAtom->GetWavePacket(otherWavepacketIndex);
        // get the position of the other wavepacket
        doubleThreeVector otherPosition = otherWavePacket->GetPosition();
        // calculate the distance between the two wavepackets
        double distance = sqrt(pow(currentPosition[0] - otherPosition[0], 2) + 
                                pow(currentPosition[1] - otherPosition[1], 2) + 
                                pow(currentPosition[2] - otherPosition[2], 2));

        // check if the distance is smaller than the coherence length
        if(distance < coherenceLength && currentWavePacket->GetState() == otherWavePacket->GetState() && wavePacketIndex != otherWavepacketIndex)
        {
            // if so, add the wavepacket indices to the group
            wavepacketIndexGroup.push_back(otherWavepacketIndex);
        }
    }

    // sort before returning
    std::sort(wavepacketIndexGroup.begin(), wavepacketIndexGroup.end());

    return wavepacketIndexGroup;
}

void AISPortFrame::initializePortVector(int numberOfPorts)
{
    fpPortVector = std::make_unique<portVector>();
    for(int portIndex = 0; portIndex < numberOfPorts; portIndex++)
    {
        std::unique_ptr<AISPort> newPort(new AISPort());
        fpPortVector->push_back(std::move(newPort));
    }

}

/*
void AISPortFrame::associatePortsWithWavePackets(AISAtom* anAtom, std::vector<std::vector<int>> wavepacketIndexGroups)
{
    int portIndex = 0;
    for(int wavePacketIndex = 0; wavePacketIndex < anAtom->GetNumberOfWavePackets() && portIndex < GetNumberOfPorts(); wavePacketIndex++)
    {
        if(isValueInVector(wavePacketIndex, adjacentWavepacketIndices[0]))
        {
            int otherWavePacketIndex = adjacentWavepacketIndices[1][findValueInVector(wavePacketIndex, adjacentWavepacketIndices[0])];
            AISPort* port = fpPortVector->at(portIndex);
            port->wavePacketIndices.push_back(wavePacketIndex);
            port->wavePacketIndices.push_back(otherWavePacketIndex);

            portIndex += 1;
        }
        else
        {
            AISPort* port = fpPortVector->at(portIndex);
            port->wavePacketIndices.push_back(wavePacketIndex);

            portIndex += 1;
        }
    
    }
*/

void AISPortFrame::setPortParameters(std::unique_ptr<AISAtom>& anAtom)
{
    // loop over ports
    for(int portIndex = 0; portIndex < GetNumberOfPorts(); portIndex++)
    {
        // get the port
        std::unique_ptr<AISPort>& port = GetPort(portIndex);
        // if only one wavepacket is in the port, no interference happens
        if(port->getNumberOfWavePackets() == 1)
        {   
            // set the interfering flag
            port->interfering = false;
            // get the wavepacket
            std::unique_ptr<AISWavePacket>& wavePacket = anAtom->GetWavePacket(port->wavePacketIndices[0]);
            // calculate the port parameters
            port->state = wavePacket->GetState();
            port->position = wavePacket->GetPosition();
            port->velocity = wavePacket->GetVelocity();

            port->phaseShift = 0.0;
            port->phaseShiftError = 0.0;

            port->probabilityAmplitude = pow(wavePacket->GetAmplitude(), 2);
        }
        // otherwise, interference happens
        else if(port->getNumberOfWavePackets() > 1)
        {   
            // set the interfering flag
            port->interfering = true;

            // check that all the wavepackets are in the same state
            int state = anAtom->GetWavePacket(port->wavePacketIndices[0])->GetState();
            for(int wavePacketIndex : port->wavePacketIndices)
            {
                if(anAtom->GetWavePacket(wavePacketIndex)->GetState() != state)
                {
                    std::cout << "Error: wavepackets in the same port are not in the same state" << std::endl;
                    exit(1);
                }
            }

            // compute the diagonal contribution to the probability amplitude
            double diagonalContribution = 0.0;
            for(int wavePacketIndex : port->wavePacketIndices)
            {
                diagonalContribution += pow(anAtom->GetWavePacket(wavePacketIndex)->GetAmplitude(), 2);
            }

            // compute the off-diagonal contribution to the probability amplitude
            double offDiagonalContribution = 0.0;
            for(int i = 0; i < port->getNumberOfWavePackets(); i++)
            {
                for(int j = i + 1; j < port->getNumberOfWavePackets(); j++)
                {
                    std::unique_ptr<AISWavePacket>& wavePacket1 = anAtom->GetWavePacket(port->wavePacketIndices[i]);
                    std::unique_ptr<AISWavePacket>& wavePacket2 = anAtom->GetWavePacket(port->wavePacketIndices[j]);

                    double A1 = wavePacket1->GetAmplitude();
                    double A2 = wavePacket2->GetAmplitude();

                    doubleThreeVector r1 = wavePacket1->GetPosition();
                    doubleThreeVector r2 = wavePacket2->GetPosition();
                    doubleThreeVector v1 = wavePacket1->GetVelocity();
                    doubleThreeVector v2 = wavePacket2->GetVelocity();

                    double phi1 = wavePacket1->GetPhaseDouble();
                    double phi2 = wavePacket2->GetPhaseDouble();

                    __float128 phi1Quad = wavePacket1->GetPhaseQuad();
                    __float128 phi2Quad = wavePacket2->GetPhaseQuad();

                    doubleThreeVector deltaX = matrixAdd(r2, scalarMultiply(r1, -1.0)); // note the inverted order
                    doubleThreeVector meanV  = matrixAdd(v1, v2);
                    meanV = scalarMultiply(meanV, 0.5);
                    double separationPhase = dotProduct(meanV, deltaX) / hbar * massSr87;
                    double dphi = phi1 - phi2 + phi1Quad - phi2Quad + separationPhase;

                    offDiagonalContribution += 2*A1*A2*cos(dphi);
                }
            }

            // compute the mean velocity and position of the group
            doubleThreeVector meanVelocity = {0.0, 0.0, 0.0};
            doubleThreeVector meanPosition = {0.0, 0.0, 0.0};
            for(int wavePacketIndex : port->wavePacketIndices)
            {
                std::unique_ptr<AISWavePacket>& wavePacket = anAtom->GetWavePacket(wavePacketIndex);
                meanVelocity = matrixAdd(meanVelocity, wavePacket->GetVelocity());
                meanPosition = matrixAdd(meanPosition, wavePacket->GetPosition());
            }
            meanVelocity = scalarMultiply(meanVelocity, 1.0/port->getNumberOfWavePackets());
            meanPosition = scalarMultiply(meanPosition, 1.0/port->getNumberOfWavePackets());

            // calculate the port parameters
            port->state = state;
            port->position = meanPosition; // take average
            port->velocity = meanVelocity; // take average
            port->probabilityAmplitude = diagonalContribution + offDiagonalContribution;

            // if the number of wavepackets is 2, the phase shift is the difference in the phase of the two wavepackets
            if(port->getNumberOfWavePackets() == 2)
            {
                std::unique_ptr<AISWavePacket>& wavePacket1 = anAtom->GetWavePacket(port->wavePacketIndices[0]);
                std::unique_ptr<AISWavePacket>& wavePacket2 = anAtom->GetWavePacket(port->wavePacketIndices[1]);
                doubleThreeVector r1 = wavePacket1->GetPosition();
                doubleThreeVector r2 = wavePacket2->GetPosition();
                doubleThreeVector v1 = wavePacket1->GetVelocity();
                doubleThreeVector v2 = wavePacket2->GetVelocity();
                double phi1 = wavePacket1->GetPhaseDouble();
                double phi2 = wavePacket2->GetPhaseDouble();
                double error1 = wavePacket1->GetPhaseDoubleError();
                double error2 = wavePacket2->GetPhaseDoubleError();
                __float128 phi1Quad = wavePacket1->GetPhaseQuad();
                __float128 phi2Quad = wavePacket2->GetPhaseQuad();

                doubleThreeVector deltaX = matrixAdd(r2, scalarMultiply(r1, -1.0)); // note the inverted order
                doubleThreeVector meanV  = matrixAdd(v1, v2);
                meanV = scalarMultiply(meanV, 0.5);
                double dphi = phi1 - phi2 + (double)(phi1Quad - phi2Quad) + dotProduct(meanV, deltaX) / hbar * massSr87;
                double dphiError = sqrt(pow(error1, 2) + pow(error2, 2));
                port->phaseShift = dphi;
                port->phaseShiftError = dphiError;

            } else
            {
                std::cout << "Error: more than two wavepackets in a port" << std::endl;
                std::cout << "Number of wavepackets: " << port->getNumberOfWavePackets() << std::endl;
                // print the median amplitude in the port
                std::vector<double> amplitudes;
                for(int wavePacketIndex : port->wavePacketIndices)
                {
                    amplitudes.push_back(anAtom->GetWavePacket(wavePacketIndex)->GetAmplitude());
                }
                std::sort(amplitudes.begin(), amplitudes.end());
                std::cout << "Median amplitude: " << amplitudes[amplitudes.size()/2] << std::endl;
            }

        } 
    }
}

int AISPortFrame::GetNumberOfPorts()
{
    return fpPortVector->size();
}

std::unique_ptr<AISPort>& AISPortFrame::GetPort(int index)
{
    return fpPortVector->at(index);
}

bool AISPortFrame::isValueInVector(int j, const intVector& vector) {
    return std::find(vector.begin(), vector.end(), j) != vector.end();
}

int AISPortFrame::findValueInVector(int j, const intVector& vector) {
    return std::find(vector.begin(), vector.end(), j) - vector.begin();
}

std::vector<std::string> AISPortFrame::GetInterferingPaths()
{
    return interferingPaths;
}
