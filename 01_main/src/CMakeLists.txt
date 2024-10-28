# Add the source files
set(SOURCES
    AISAtom.cc
    AISAtomEnsemble.cc
    AISDataIO.cc
    AISDetector.cc
    AISDriver.cc
    AISKinematicPropagator.cc
    AISLaserBeam.cc
    AISParams.cc
    AISPort.cc
    AISPortFrame.cc
    AISPulsePropagator.cc
    AISUtilities.cc
    AISWavePacket.cc
)

# Add the library target
add_library(AISppLib ${SOURCES})

# Link the library with the required libraries
target_link_libraries(AISppLib AISPotentials AISBeams ${HDF5_LIBRARIES} ${GSL_LIBRARIES} OpenMP::OpenMP_CXX quadmath)

# Include directories for the library
target_include_directories(AISppLib PUBLIC ${CMAKE_SOURCE_DIR}/01_main/include)