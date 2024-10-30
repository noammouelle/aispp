cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

project(AISppLib VERSION 0.1 LANGUAGES CXX)

# Set the build type to Debug if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

# Find the required libraries
link_directories(/lib64/libstdc++.so.6)
#set(HDF5_ROOT "/work/ndm33/local/hdf5")
set(HDF5_INCLUDE_DIRS "${HDF5_ROOT}/include")
set(HDF5_LIBRARIES "${HDF5_ROOT}/lib/libhdf5.so")
find_package(HDF5 REQUIRED COMPONENTS CXX)
find_package(OpenMP REQUIRED)
find_package(GSL REQUIRED)

# Set the include directories
include_directories(${HDF5_INCLUDE_DIRS}) # HDF5 libraries
include_directories(${GSL_INCLUDE_DIRS})  # GSL libraries
include_directories(${CMAKE_SOURCE_DIR}/01_main/include) # Project headers

# Add OpenMP flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -O0")

# Add the source directories
add_subdirectory(01_main/src)
add_subdirectory(02_aispotentials)
add_subdirectory(03_aisbeams)

# Specify the executable target
add_executable(ais++ 01_main/src/AIS++.cc)

# Link the executable with the required libraries
target_link_libraries(ais++ AISppLib AISPotentials AISBeams ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${GSL_LIBRARIES} OpenMP::OpenMP_CXX quadmath)