# AIS++ v2.0.2

Atom Interferometry Simulator in C++ (AIS++) is a C++ library enabling the simulation of atom clouds undergoing atom interferometry sequences. 

AIS++ version > 1.0.0 uses numerical ODE solvers to solve for the kinematics of the wavepackets in-between pulses, as well as for the evolution of the wavepackets during pulses. v2.0.2 implements a corrected version of the beam-splitting scheme, documented in the LaTeX proejct "Numerical treatment semi-classical beam-splitter" as well as my notebook. This version also includes several tests to ensure the correct implementation of the beam-splitter. 

## Table of Contents

- [Introduction](#introduction)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [References](#references)
- [License](#license)

## Introduction

## Dependencies
For ais++
- GNU Scientific Library 2.5 (instructions at https://coral.ise.lehigh.edu/jild13/2016/07/11/hello/)
- HDF5
- CMake
- G++
- GCC

For the examples
- Matplotlib
- h5py
- Numpy
- Pandas

## Installation
Start by cloning this repository
```
git clone git@github.com:noammouelle/aispp.git
```
You then need edit the `config.cmake` file, uncommenting the relevant lines and adding the path to your GSL and HDF5 installations. Once done, create a build directory
```
mkdir build
cd build
```
and create the make files using cmake, passing as an argument the path to you configuration file
```
cmake -DCMAKE_CONFIG_FILE=/path/to/config.cmake ..
```
Once done, build the executables using
```
make
```
Finally, you will need to add the following lines to your `.bashrc` file
```
export AISPP_BUILD="/path/to/ais++/build"
export PATH="$AISPP_BUILD:$PATH"
```

## Usage
The code is run by calling the executable `ais++` from terminal, specifying the input and output files as such
```
ais++ -i /path/to/input.aisi -o /path/to/output.h5
```


## Contributing

Email: ndm33@cam.ac.uk

## References


## Licence
