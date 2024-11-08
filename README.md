# AIS++ v2.0.0

Atom Interferometry Simulator in C++ (AIS++) is a C++ library enabling the simulation of atom clouds undergoing atom interferometry sequences. 

In v2.0.0, a new scheme is used to model the propagation of the atoms in laser beams. Also, the code relies on numerical ODE solvers and numerical integration to propagate the atoms and compute classical actions along the atoms trajectories. This version also intends to be more modular, general, and to make calculations more tractable. Potentials, and sequence parameters are now defined externally by the user in a .aisi configuration file, avoiding the need to recompile the code for different configurations.

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
- GNU Scientific Library 2.5 (instructions at https://coral.ise.lehigh.edu/jild13/2016/07/11/hello/)
- HDF5


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
