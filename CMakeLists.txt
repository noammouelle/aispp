cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

project(AISppLib VERSION 0.1 LANGUAGES CXX)

# Set the build type to Debug if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

# Find the required libraries
link_directories(/lib64/libstdc++.so.6)
set(HDF5_ROOT "/work/ndm33/local/hdf5")
#set(HDF5_INCLUDE_DIRS "${HDF5_ROOT}/include")
#set(HDF5_LIBRARIES "${HDF5_ROOT}/lib/libhdf5.so")
find_package(HDF5 REQUIRED COMPONENTS CXX)
find_package(OpenMP REQUIRED)
find_package(GSL REQUIRED)


# Set the include directories
include_directories(${HDF5_INCLUDE_DIRS}) # HDF5 libraries
include_directories(${GSL_INCLUDE_DIRS})  # GSL libraries

# Add OpenMP flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# Set compiler flags
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    # Debugging flags: disable optimizations and enable debug symbols
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g")
else()
    # Release mode: add optimization flags
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")
endif()

# Collect the source files
file(GLOB SOURCES ${PROJECT_SOURCE_DIR}/src/*.cc) # all source files
file(GLOB HEADERS ${PROJECT_SOURCE_DIR}/include/*.hh) # all header files

# Create a library
add_library(AISppLib ${SOURCES})

# Link the libraries
target_link_libraries(AISppLib PRIVATE ${HDF5_LIBRARIES} OpenMP::OpenMP_CXX ${GSL_LIBRARIES})

# Set the include directories for the library
target_include_directories(AISppLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Optional: Define GSL_DEBUG for debugging mode
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(AISppLib PRIVATE GSL_DEBUG)
endif()
