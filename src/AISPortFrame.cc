#include "AISPortFrame.hh"

AISPortFrame::AISPortFrame(AISAtom* anAtom, intTuple adjacentWavepacketIndices)
{
    // create a flag vector indicating if a wavepacket must be accounted for or not
    boolVector considerWavepacket(anAtom->GetNumberOfWavePackets(), true);
    // loop over all the wavepackets of the atom
    for(int wavepacketIndex = 0; wavepacketIndex < anAtom->GetNumberOfWavePackets(); wavepacketIndex++)
    {
        // check if it should be considered
        if(considerWavepacket[wavepacketIndex])
        {
            // create a new port
            AISPort* newPort = new AISPort();
            // check if the wavepacket is in the list of adjacent wavepackets
            if(isValueInVector(wavepacketIndex, adjacentWavepacketIndices[0]) == false && isValueInVector(wavepacketIndex, adjacentWavepacketIndices[1]) == false)
            {
                // if not, the port properties are defined as such
                AISWavePacket* currentWavePacket = anAtom->GetWavePacket(wavepacketIndex);
                newPort->state = currentWavePacket->GetState();
                newPort->probabilityAmplitude = pow(currentWavePacket->GetAmplitude(), 2);
                newPort->position = currentWavePacket->GetPosition();
                newPort->velocity = currentWavePacket->GetVelocity();
                newPort->phaseShift = 0.0;
                // add the port to the port vector
                fpPortVector->push_back(newPort);
            } else if(isValueInVector(wavepacketIndex, adjacentWavepacketIndices[0]) == true)
            {
                // find the index of the wavepacket index in the adjacent wavepacket indices list
                int index = std::find(adjacentWavepacketIndices[0].begin(), adjacentWavepacketIndices[0].end(), wavepacketIndex) - adjacentWavepacketIndices[0].begin();
                AISWavePacket* wavePacket1 = anAtom->GetWavePacket(wavepacketIndex);
                AISWavePacket* wavePacket2 = anAtom->GetWavePacket(adjacentWavepacketIndices[1][index]);

                newPort->state = wavePacket1->GetState();
                newPort->position = scalarMultiply((wavePacket1->GetPosition(), wavePacket2->GetPosition()), 0.5);
                newPort->velocity = scalarMultiply((wavePacket1->GetVelocity(), wavePacket2->GetVelocity()), 0.5);
                newPort->phaseShift = wavePacket1->GetPhaseDouble() - wavePacket2->GetPhaseDouble() 
                                      + convertScalarToDouble(wavePacket1->GetPhaseQuad() - wavePacket2->GetPhaseQuad());

                doubleThreeVector deltaX = matrixAdd(wavePacket2->GetPosition(), 
                                                     scalarMultiply(wavePacket1->GetPosition(), -1.0)); // note the inverted order
                
                newPort->phaseShift += dotProduct(newPort->velocity, deltaX) / hbar * massSr87;

                // account for interference in probability calculation
                double A1 = wavePacket1->GetAmplitude();
                double A2 = wavePacket2->GetAmplitude();
                double dphi = newPort->phaseShift;
                newPort->probabilityAmplitude = pow(A1,2) + pow(A2,2) + 2*A1*A2*cos(dphi);
                
                // add the port to the port vector
                fpPortVector->push_back(newPort);

                // mark the wavepackets as considered
                considerWavepacket[wavepacketIndex] = false;
                considerWavepacket[adjacentWavepacketIndices[1][index]] = false;

            }
        }
    }
}

AISPortFrame::~AISPortFrame()
{
    for(AISPort* port : *fpPortVector)
    {
        delete port;
    }
    delete fpPortVector;
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