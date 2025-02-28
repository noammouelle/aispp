cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

project(AISppLib VERSION 0.1 LANGUAGES CXX)

# Check if a configuration file is provided
if(NOT DEFINED CMAKE_CONFIG_FILE)
    message(FATAL_ERROR "CMAKE_CONFIG_FILE is not defined. Please specify the configuration file using -DCMAKE_CONFIG_FILE=<path_to_config_file>")
endif()

# Include the configuration file
include(${CMAKE_CONFIG_FILE})

# Find the required libraries
find_package(HDF5 REQUIRED COMPONENTS CXX)
find_package(OpenMP REQUIRED)
find_package(GSL REQUIRED)

# Set the include directories
include_directories(${HDF5_INCLUDE_DIRS}) # HDF5 libraries
include_directories(${GSL_INCLUDE_DIRS})  # GSL libraries
include_directories(${EIGEN3_INCLUDE_DIR}) # Eigen library
include_directories(${CMAKE_SOURCE_DIR}/01_main/include) # Project headers

# Add OpenMP flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -O3 -march=native -pg")

# Add the source directories
add_subdirectory(01_main/src)
add_subdirectory(02_aispotentials)
add_subdirectory(03_aisbeams)

# Specify the executable target
add_executable(ais++ 01_main/src/AIS++.cc)

# Link the executable with the required libraries
target_link_libraries(ais++ AISppLib AISPotentials AISBeams ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${GSL_LIBRARIES} OpenMP::OpenMP_CXX quadmath)