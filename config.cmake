# Template configuration file for ais++

# path to the HDF5 library (optional)
# set(HDF5_ROOT "/work/ndm33/local/hdf5")

# path to the GSL library (optional)
# set(GSL_ROOT "/work/ndm33/local/gsl")

# Set the build type to Debug if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()