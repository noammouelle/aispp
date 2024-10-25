#ifndef AISCONSTANTS_HH
#define AISCONSTANTS_HH

// Defines the physical constants required for the
// simulation.

#include <cmath>

const double g    = 9.81;
const double h    = 6.62607015e-34;
const double pi   = 3.141592653589793;
const double hbar = h/(2*pi);
const double kB   = 1.3806488e-23;
const int   c     = 299792458;

const __float128 pi128   = 3.1415926535897932384626433832795028q;

const __float128 frequencySr87 = 429228004229873.0q;
const __float128 omegaSr87 = 2.0q * pi128 * frequencySr87;
const __float128 lambdaSr87_quad = c * 1.0q / frequencySr87;

const double au          = 1.660539066e-27;
const double massSr87    = 86.90888 * au;
const double lambdaSr87  = static_cast<double>(lambdaSr87_quad);

const double massEarth   = 5.97e24;
const double radiusEarth = 6.37e6;

#endif