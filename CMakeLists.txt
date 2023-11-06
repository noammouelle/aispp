cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

project(AISim++)

# Find the HDF5 library
find_package(HDF5 REQUIRED COMPONENTS CXX)
include_directories(${HDF5_INCLUDE_DIRS})

include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/src)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc) # all source files
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(aispp src/main.cc ${sources} ${headers})

target_link_libraries(aispp ${HDF5_LIBRARIES})

add_custom_target(AISim++ DEPENDS aispp)