cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

project(AISppLib VERSION 0.1 LANGUAGES CXX)

# Find the required libraries
link_directories(/lib64/libstdc++.so.6)
set(HDF5_ROOT "/work/ndm33/local/hdf5")
#set(HDF5_INCLUDE_DIRS "${HDF5_ROOT}/include")
#set(HDF5_LIBRARIES "${HDF5_ROOT}/lib/libhdf5.so")
find_package(HDF5 REQUIRED COMPONENTS CXX)
find_package(OpenMP REQUIRED)


# Set the include directories
include_directories(${HDF5_INCLUDE_DIRS}) # libraries

#include_directories(${CMAKE_SOURCE_DIR}/include)
#include_directories(${CMAKE_SOURCE_DIR}/src)

# Add OpenMP flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# Set compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

# Collect the source files
file(GLOB SOURCES ${PROJECT_SOURCE_DIR}/src/*.cc) # all source files
#file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh) # all header files

# create a library
add_library(AISppLib ${SOURCES})

# Link the libraries
target_link_libraries(AISppLib PRIVATE ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES}  OpenMP::OpenMP_CXX)

# Set the include directories
target_include_directories(AISppLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
