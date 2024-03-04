#include "AISPortFrame.hh"

AISPortFrame::AISPortFrame(AISAtom* anAtom, intTuple adjacentWavepacketIndices)
{
    // count how many ports are needed
    int numberOfPorts = anAtom->GetNumberOfWavePackets() - adjacentWavepacketIndices[0].size();

    // error if the number of ports is non-positive
    if(numberOfPorts <= 0)
    {
        std::cout << "ERROR: number of ports must be positive" << std::endl;
        std::cout << "Number of ports: " << numberOfPorts << std::endl;
        std::cout << "Number of wavepackets: " << anAtom->GetNumberOfWavePackets() << std::endl;
        std::cout << "Number of adjacent wavepackets: " << adjacentWavepacketIndices[0].size() << std::endl;
        exit(1);
    }

    // create the port vector
    initializePortVector(numberOfPorts);
    // associate the ports with the wavepackets
    associatePortsWithWavePackets(anAtom, adjacentWavepacketIndices);
    // set the port parameters
    setPortParameters(anAtom, adjacentWavepacketIndices);
}

AISPortFrame::~AISPortFrame()
{
    for(AISPort* port : *fpPortVector)
    {
        delete port;
    }
    delete fpPortVector;
}

void AISPortFrame::initializePortVector(int numberOfPorts)
{
    fpPortVector = new portVector();
    for(int portIndex = 0; portIndex < numberOfPorts; portIndex++)
    {
        AISPort* newPort = new AISPort();
        fpPortVector->push_back(newPort);
    }

}

void AISPortFrame::associatePortsWithWavePackets(AISAtom* anAtom, intTuple adjacentWavepacketIndices)
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
}

void AISPortFrame::setPortParameters(AISAtom* anAtom, intTuple adjacentWavepacketIndidces)
{
    // loop over ports
    for(int portIndex = 0; portIndex < GetNumberOfPorts(); portIndex++)
    {
        // get the port
        AISPort* port = GetPort(portIndex);
        // check if two wavepackets interfere
        if(port->getNumberOfWavePackets() == 2)
        {
            // get the wavepackets
            AISWavePacket* wavePacket1 = anAtom->GetWavePacket(port->wavePacketIndices[0]);
            AISWavePacket* wavePacket2 = anAtom->GetWavePacket(port->wavePacketIndices[1]);

            int state1 = wavePacket1->GetState();
            int state2 = wavePacket2->GetState();
            assert(state1 == state2);

            doubleThreeVector r1 = wavePacket1->GetPosition();
            doubleThreeVector r2 = wavePacket2->GetPosition();
            doubleThreeVector v1 = wavePacket1->GetVelocity();
            doubleThreeVector v2 = wavePacket2->GetVelocity();

            double phi1 = wavePacket1->GetPhaseDouble();
            double phi2 = wavePacket2->GetPhaseDouble();

            __float128 phi1Quad = wavePacket1->GetPhaseQuad();
            __float128 phi2Quad = wavePacket2->GetPhaseQuad();

            // calculate the port parameters
            port->state = state1;
            port->position = scalarMultiply(matrixAdd(r1, r2), 0.5); // take average
            port->velocity = scalarMultiply(matrixAdd(v1, v2), 0.5);
            port->phaseShift = phi1 - phi2 + phi1Quad - phi2Quad;

            doubleThreeVector deltaX = matrixAdd(r2, scalarMultiply(r1, -1.0)); // note the inverted order
            port->phaseShift += dotProduct(port->velocity, deltaX) / hbar * massSr87;

            // account for interference in probability calculation
            double A1 = wavePacket1->GetAmplitude();
            double A2 = wavePacket2->GetAmplitude();
            double dphi = port->phaseShift;
            port->probabilityAmplitude = pow(A1,2) + pow(A2,2) + 2*A1*A2*cos(dphi);

        } else if(port->getNumberOfWavePackets() == 1)
        {
            // get the wavepacket
            AISWavePacket* wavePacket = anAtom->GetWavePacket(port->wavePacketIndices[0]);
            // calculate the port parameters
            port->state = wavePacket->GetState();
            port->position = wavePacket->GetPosition();
            port->velocity = wavePacket->GetVelocity();

            port->phaseShift = 0.0;

            port->probabilityAmplitude = pow(wavePacket->GetAmplitude(), 2);
        }
    }
}

int AISPortFrame::GetNumberOfPorts()
{
    return fpPortVector->size();
}

AISPort* AISPortFrame::GetPort(int index)
{
    return fpPortVector->at(index);
}

bool AISPortFrame::isValueInVector(int j, const intVector& vector) {
    return std::find(vector.begin(), vector.end(), j) != vector.end();
}

int AISPortFrame::findValueInVector(int j, const intVector& vector) {
    return std::find(vector.begin(), vector.end(), j) - vector.begin();
}