#ifndef AISDATAIO_HH
#define AISDATAIO_HH

#include <iostream>
#include <array>
#include <vector>
#include <cmath>
#include <string>
#include <complex>
#include <sstream>
//#include <quadmath.h>
#include <iomanip>
#include <iomanip>
#include <fstream>

#include "AISConstants.hh"
#include "AISAtomEnsemble.hh"

/* I/O routines */
void writeAtomEnsembleToFile(std::string fName, AISAtomEnsemble* atomEnsemble);


#endif