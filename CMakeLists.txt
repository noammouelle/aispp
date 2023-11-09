cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

project(AISim++)

# Find the required libraries
find_package(HDF5 REQUIRED COMPONENTS CXX)
find_package(OpenMP REQUIRED)

# Set the include directories
include_directories(${HDF5_INCLUDE_DIRS}) # libraries

include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/src)

# Add OpenMP flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc) # all source files
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(aispp src/main.cc ${sources} ${headers})

target_link_libraries(aispp ${HDF5_LIBRARIES} OpenMP::OpenMP_CXX)

add_custom_target(AISim++ DEPENDS aispp)