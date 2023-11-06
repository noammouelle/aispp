#include "AISDetector.hh"

AISDetector::AISDetector(AISAtomEnsemble* pAtomEnsemble, double coherenceLength)
{
    for(int atomIndex = 0; atomIndex < pAtomEnsemble->GetNumberOfAtoms(); atomIndex++)
    {
        AISAtom* currentAtom = pAtomEnsemble->GetAtom(atomIndex);
        // get the indices of waveapckets close enough to interfere
        intTuple adjacentWavepacketIndices = GetAdjacentWavepackets(currentAtom, coherenceLength);
        // create the port-frame object
        AISPortFrame* currentPortFrame = new AISPortFrame(currentAtom, adjacentWavepacketIndices);
        fpPortFrameVector->push_back(currentPortFrame);
    }
}

AISDetector::~AISDetector()
{
    delete fpPortFrameVector;
}

intTuple AISDetector::GetAdjacentWavepackets(AISAtom* anAtom, double coherenceLength)
{
    intTuple adjacentWavepacketIndices;
    // loop over all the wavepackets of the atom
    for(int wavepacketIndex = 0; wavepacketIndex < anAtom->GetNumberOfWavePackets(); wavepacketIndex++)
    {
        // get the current wavepacket
        AISWavePacket* currentWavePacket = anAtom->GetWavePacket(wavepacketIndex);
        // get the position of the current wavepacket
        doubleThreeVector currentPosition = currentWavePacket->GetPosition();
        // loop over all the other wavepackets
        for(int otherWavepacketIndex = wavepacketIndex + 1; otherWavepacketIndex < anAtom->GetNumberOfWavePackets(); otherWavepacketIndex++)
        {
            // get the other wavepacket
            AISWavePacket* otherWavePacket = anAtom->GetWavePacket(otherWavepacketIndex);
            // get the position of the other wavepacket
            doubleThreeVector otherPosition = otherWavePacket->GetPosition();
            // calculate the distance between the two wavepackets
            double distance = sqrt(pow(currentPosition[0] - otherPosition[0], 2) + 
                                   pow(currentPosition[1] - otherPosition[1], 2) + 
                                   pow(currentPosition[2] - otherPosition[2], 2));
            // check if the distance is smaller than the coherence length
            if(distance < coherenceLength)
            {
                // if so, add the wavepacket indices to the list
                adjacentWavepacketIndices[0].push_back(wavepacketIndex);
                adjacentWavepacketIndices[1].push_back(otherWavepacketIndex);
            }
        }
    }
    return adjacentWavepacketIndices;
}

int AISDetector::SamplePort(AISPortFrame* aPortFrame)
{
    // compute the cumulative probabilities
    doubleVector cumulativeProbabilities;
    double cumulativeProbability = 0.0;
    for(int portIndex = 0; portIndex < aPortFrame->GetNumberOfPorts(); portIndex++)
    {
        cumulativeProbability += aPortFrame->GetPort(portIndex)->probabilityAmplitude;
        cumulativeProbabilities.push_back(cumulativeProbability);
    }

    // sample a random number
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    double randomNumber = dis(gen);

    // find the port corresponding to the random number
    // Find the index of the sampled object
    for (int i = 0; i < cumulativeProbabilities.size(); ++i) {
        if (randomNumber <= cumulativeProbabilities[i]) {
            return i;
        }
    }

    // Handle an error case where no object is sampled
    return -1;
}

intVector AISDetector::SampleAllPorts()
{
    intVector sampledPortsIndices;
    // loop over all the port-frames
    numSamples = 0;
    for(int portFrameIndex = 0; portFrameIndex < GetNumberOfPortFrames(); portFrameIndex++)
    {
        // sample a port from the current port-frame
        int sampledPortIndex = SamplePort(GetPortFrame(portFrameIndex));
        // add the sampled port index to the list
        sampledPortsIndices.push_back(sampledPortIndex);
        // increment the number of samples if index is not -1
        if(sampledPortIndex != -1)
        {
            numSamples++;
        }
    }
    return sampledPortsIndices;
}

int AISDetector::GetNumberOfPortFrames()
{
    return fpPortFrameVector->size();
}

int AISDetector::GetNumberOfSamples()
{
    return numSamples;
}

AISPortFrame* AISDetector::GetPortFrame(int portFrameIndex)
{
    return fpPortFrameVector->at(portFrameIndex);
}