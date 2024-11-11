#include "AISDriver.hh"
#include "AISDataIO.hh"
#include "AISParams.hh"

#include <iostream>
#include <chrono>

/*
Main file for this library. Works by calling 
    ais++ -i input_file.aisi -o output_file.h5
in the terminal. The input file should be a .aisi file
*/

int main(int argc, char* argv[])
{
    // parse the command line arguments
    std::string inputFileName;
    std::string outputFileName;
    if(argc == 5)
    {
        if(std::string(argv[1]) == "-i" && std::string(argv[3]) == "-o")
        {
            inputFileName = std::string(argv[2]);
            outputFileName = std::string(argv[4]);
        }
        else
        {
            std::cerr << "Usage: ais++ -i input_file.aisi -o output_file.h5" << std::endl;
            exit(1);
        }
    }
    else
    {
        std::cerr << "Usage: ais++ -i input_file.aisi -o output_file.h5" << std::endl;
        exit(1);
    }

    // read the params
    AISParams params = readParamsFromFile(inputFileName);
    
    // create the driver
    AISDriver* driver = new AISDriver(params);
    
    /* RUN THE ATOM PROPAGATION */
    auto start = std::chrono::high_resolution_clock::now();

    driver->Run();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Propagation: " << elapsed.count() << " s" << std::endl;

    /* DETECT THE ATOMS */
    auto startDetect = std::chrono::high_resolution_clock::now();

    driver->Detect();

    auto endDetect = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedDetect = endDetect - startDetect;
    std::cout << "Detection: " << elapsedDetect.count() << " s" << std::endl;

    // write the detected atoms to file
    auto startWrite = std::chrono::high_resolution_clock::now();

    driver->WriteDetectedAtomsToFile(outputFileName);

    // if printPorts = true, print the ports too
    if (driver->params.printPorts) {
    std::string portFileName = outputFileName.substr(0, outputFileName.find_last_of('.')) + "_PROB.h5";
    driver->WritePortsToFile(portFileName);
    }

    auto endWrite = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedWrite = endWrite - startWrite;
    std::cout << "Write: " << elapsedWrite.count() << " s" << std::endl;

    // clean up
    delete driver;

    return 0;
}