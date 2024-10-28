#include "AISDriver.hh"
#include "AISDataIO.hh"
#include "AISParams.hh"

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

    // run the simulation
    driver->Run();
    driver->Detect();

    // write the detected atoms to file
    driver->WriteDetectedAtomsToFile(outputFileName);

    // clean up
    delete driver;

    return 0;
}