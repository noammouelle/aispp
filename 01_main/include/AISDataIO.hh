#ifndef AISDATAIO_HH
#define AISDATAIO_HH

#include <iostream>
#include <fstream>
#include <sstream>
#include <array>
#include <vector>
#include <set>
#include <map>
#include <cmath>
#include <string>
#include <complex>
#include <sstream>
#include <quadmath.h>
#include <iomanip>
#include <iomanip>
#include <fstream>

#include "AISConstants.hh"
#include "AISAtomEnsemble.hh"
#include "AISParams.hh"

/* I/O routines */
void writeAtomEnsembleToFile(std::string fName, std::unique_ptr<AISAtomEnsemble>& atomEnsemble);
AISParams readParamsFromFile(std::string fName);


#endif