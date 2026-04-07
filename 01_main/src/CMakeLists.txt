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

# GPU source — only compiled when USE_CUDA is ON
if(USE_CUDA)
    list(APPEND SOURCES AISPulsePropagatorGPU.cu AISKinematicPropagatorGPU.cu)
endif()

# Add the library target
add_library(AISppLib ${SOURCES})

# CUDA requires separable compilation when device code is split across files
if(USE_CUDA)
    set_target_properties(AISppLib PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    target_link_libraries(AISppLib AISPotentials AISBeams ${HDF5_LIBRARIES} ${GSL_LIBRARIES} OpenMP::OpenMP_CXX quadmath CUDA::cudart)
else()
    target_link_libraries(AISppLib AISPotentials AISBeams ${HDF5_LIBRARIES} ${GSL_LIBRARIES} OpenMP::OpenMP_CXX quadmath)
endif()

# Include directories for the library
target_include_directories(AISppLib PUBLIC ${CMAKE_SOURCE_DIR}/01_main/include)