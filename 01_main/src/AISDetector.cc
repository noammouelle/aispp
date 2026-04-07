#include "AISDetector.hh"

AISDetector::AISDetector(std::unique_ptr<AISAtomEnsemble>& pAtomEnsemble, double coherenceLength)
{
    int nAtoms = pAtomEnsemble->GetNumberOfAtoms();

    // Pre-allocate the vector so each parallel thread can write to its own
    // index without any synchronisation.  The previous push_back inside an
    // omp critical section serialised all 1M insertions, making detection
    // effectively single-threaded despite the parallel for.
    fpPortFrameVector = std::make_unique<portFrameVector>(nAtoms);

    #pragma omp parallel for
    for(int atomIndex = 0; atomIndex < nAtoms; atomIndex++)
    {
        std::unique_ptr<AISAtom>& currentAtom = pAtomEnsemble->GetAtom(atomIndex);
        (*fpPortFrameVector)[atomIndex] =
            std::make_unique<AISPortFrame>(currentAtom, coherenceLength);
    }
}

AISDetector::~AISDetector()
{
    // delete all the unique pointers to port-frames then free the vector smart pointer
    fpPortFrameVector->clear();
    fpPortFrameVector.reset();
}

/*
intTuple AISDetector::GetAdjacentWavepackets(std::unique_ptr<AISAtom>& anAtom, double coherenceLength)
{
    intTuple adjacentWavepacketIndices;
    // loop over all the wavepackets of the atom
    for(int wavepacketIndex = 0; wavepacketIndex < anAtom->GetNumberOfWavePackets(); wavepacketIndex++)
    {
        // get the current wavepacket
        std::unique_ptr<AISWavePacket>& currentWavePacket = anAtom->GetWavePacket(wavepacketIndex);
        // get the position of the current wavepacket
        doubleThreeVector currentPosition = currentWavePacket->GetPosition();
        // loop over all the other wavepackets
        for(int otherWavepacketIndex = wavepacketIndex + 1; otherWavepacketIndex < anAtom->GetNumberOfWavePackets(); otherWavepacketIndex++)
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
            if(distance < coherenceLength && currentWavePacket->GetState() == otherWavePacket->GetState())
            {
                // if so, add the wavepacket indices to the list
                adjacentWavepacketIndices[0].push_back(wavepacketIndex);
                adjacentWavepacketIndices[1].push_back(otherWavepacketIndex);
            }
        }
    }
    return adjacentWavepacketIndices;
}
*/
int AISDetector::SamplePort(std::unique_ptr<AISPortFrame>& aPortFrame)
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
    sampledPortsIndices.reserve(GetNumberOfPortFrames());
    numSamples = 0;

    // The original SamplePort() created a new std::random_device and std::mt19937
    // inside the per-atom loop, causing 1M syscall-based entropy reads from
    // /dev/urandom and 1M Mersenne Twister initializations.  Seed one RNG here
    // and reuse it for all atoms.  Also avoid the per-atom cumulativeProbabilities
    // vector allocation by inlining the walk.
    std::random_device seedSource;
    std::mt19937 rng(seedSource());
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    for(int portFrameIndex = 0; portFrameIndex < GetNumberOfPortFrames(); portFrameIndex++)
    {
        auto& portFrame = GetPortFrame(portFrameIndex);
        double r = dis(rng);
        double cumulative = 0.0;
        int sampledPortIndex = -1;
        for(int portIndex = 0; portIndex < portFrame->GetNumberOfPorts(); portIndex++)
        {
            cumulative += portFrame->GetPort(portIndex)->probabilityAmplitude;
            if(r <= cumulative) { sampledPortIndex = portIndex; break; }
        }
        sampledPortsIndices.push_back(sampledPortIndex);
        if(sampledPortIndex != -1) numSamples++;
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

std::unique_ptr<AISPortFrame>& AISDetector::GetPortFrame(int portFrameIndex)
{
    // assert that the port-frame vector is not empty
    assert(fpPortFrameVector->size() > 0);
    // assert that the port-frame index is valid
    assert(portFrameIndex < GetNumberOfPortFrames());
    // return the port-frame
    return fpPortFrameVector->at(portFrameIndex);
}