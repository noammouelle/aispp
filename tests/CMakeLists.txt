# Unit tests. Built as a separate target so the normal `ais++` build is unaffected:
#   cmake --build build --target test_rotating && ./build/tests/test_rotating
add_executable(test_rotating EXCLUDE_FROM_ALL test_rotating.cc)

target_link_libraries(test_rotating
    AISppLib
    AISPotentials
    AISBeams
    ${HDF5_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
    ${GSL_LIBRARIES}
    OpenMP::OpenMP_CXX
    quadmath)

add_custom_target(check
    COMMAND test_rotating
    DEPENDS test_rotating
    COMMENT "Running ais++ unit tests")
