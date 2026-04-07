cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

# ---------------------------------------------------------------------------
# CUDA support (optional) — must be declared before the first project() call
# so CMake can detect the CUDA compiler.
# Usage: cmake -DUSE_CUDA=ON ..
# ---------------------------------------------------------------------------
option(USE_CUDA "Enable GPU acceleration via CUDA (ultrafast mode only)" OFF)

if(USE_CUDA)
    project(AISppLib VERSION 0.1 LANGUAGES CXX CUDA)
else()
    project(AISppLib VERSION 0.1 LANGUAGES CXX)
endif()

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
if(USE_CUDA)
    find_package(CUDAToolkit REQUIRED)
endif()

# Set the include directories
include_directories(${HDF5_INCLUDE_DIRS}) # HDF5 libraries
include_directories(${GSL_INCLUDE_DIRS})  # GSL libraries
include_directories(${EIGEN3_INCLUDE_DIR}) # Eigen library
include_directories(${CMAKE_SOURCE_DIR}/01_main/include) # Project headers

# Add OpenMP flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -O3 -march=native")

# CUDA-specific settings
if(USE_CUDA)
    # Propagate the USE_CUDA preprocessor define to all targets
    add_compile_definitions(USE_CUDA)
    # Default GPU architecture — override with -DCUDA_ARCH=89 for Ada, etc.
    if(NOT DEFINED CUDA_ARCH)
        set(CUDA_ARCH "75")   # Turing (RTX 20xx) — safe default for modern GPUs
    endif()
    set(CMAKE_CUDA_ARCHITECTURES ${CUDA_ARCH})
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3 --use_fast_math")
    message(STATUS "CUDA enabled — targeting sm_${CUDA_ARCH} (override with -DCUDA_ARCH=XX)")
endif()

# Add the source directories
add_subdirectory(01_main/src)
add_subdirectory(02_aispotentials)
add_subdirectory(03_aisbeams)

# Specify the executable target
add_executable(ais++ 01_main/src/AIS++.cc)

# Link the executable with the required libraries
target_link_libraries(ais++ AISppLib AISPotentials AISBeams ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${GSL_LIBRARIES} OpenMP::OpenMP_CXX quadmath)