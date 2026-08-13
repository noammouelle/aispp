# Unit tests. Built as separate targets so the normal `ais++` build is
# unaffected:
#   cmake --build build --target check
# or individually:
#   cmake --build build --target test_rotating && ./build/tests/test_rotating

set(AISPP_TEST_LIBS
    AISppLib
    AISPotentials
    AISBeams
    ${HDF5_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
    ${GSL_LIBRARIES}
    OpenMP::OpenMP_CXX
    quadmath)

# Rotating-frame (Coriolis) potentials and the velocity-dependent terms.
add_executable(test_rotating EXCLUDE_FROM_ALL test_rotating.cc)
target_link_libraries(test_rotating ${AISPP_TEST_LIBS})

# Sampled-beam path: wtype "interpolated", the aisoptics HDF5 layout, the
# trilinear interpolant, and the analytic beam types it must not disturb.
add_executable(test_interpolation EXCLUDE_FROM_ALL test_interpolation.cc)
target_link_libraries(test_interpolation ${AISPP_TEST_LIBS})

add_custom_target(check
    COMMAND test_rotating
    COMMAND test_interpolation
    DEPENDS test_rotating test_interpolation
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Running ais++ unit tests")

# The Coriolis end-to-end suite drives real Mach-Zehnder sequences through the
# built executable and needs aispy on PYTHONPATH, so it is a separate target
# rather than part of `check`:
#   cmake --build build --target check-e2e
add_custom_target(check-e2e
    COMMAND ${CMAKE_COMMAND} -E env AISPP_BINARY=$<TARGET_FILE:ais++>
            python3 ${CMAKE_CURRENT_SOURCE_DIR}/test_coriolis_e2e.py
    DEPENDS ais++
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Running ais++ Coriolis end-to-end tests (needs aispy)")
