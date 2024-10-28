# Add the source files
set(AISBEAMS_SOURCES
    src/AISEnvelopes.cc
    src/AISWavefronts
)

# Add the library target
add_library(AISBeams ${AISBEAMS_SOURCES})

# Link the library with the required libraries
target_link_libraries(AISppLib ${HDF5_LIBRARIES} ${GSL_LIBRARIES} OpenMP::OpenMP_CXX)

# Include directories for the library
target_include_directories(AISBeams PUBLIC ${CMAKE_SOURCE_DIR}/03_aisbeams/include)