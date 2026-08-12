#include "AISPotentials.hh"

// no gravity
// all in .hh now for inlining.

// uniform gravity
// all in .hh now for inlining.

// linear gravity
// all in .hh now for inlining.

// rotating frame
// The potential functions are plain function pointers, so the frame rotation
// rate cannot be captured and lives here instead. Zero means "not rotating",
// in which case every rotating_* potential reduces exactly to its inertial
// counterpart.
doubleThreeVector gRotationRate = {0.0, 0.0, 0.0};

void SetRotationRate(const doubleThreeVector& omega)
{
    gRotationRate = omega;
}
