# Add the source files
set(AISPOTENTIALS_SOURCES
    src/AISPotentials.cc
)

# Add the library target
add_library(AISPotentials ${AISPOTENTIALS_SOURCES})

# Include directories for the library
target_include_directories(AISPotentials PUBLIC ${CMAKE_SOURCE_DIR}/02_aispotentials/include)