# Copy test data to build directory
file(COPY data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Fetch googletest 
FetchContent_MakeAvailable(googletest)
set_property(DIRECTORY "${googletest_SOURCE_DIR}" PROPERTY EXCLUDE_FROM_ALL YES)

include(GoogleTest)
include(cmake/KaTestrophe.cmake)

function (kagen_add_test KAGEN_TARGET_NAME)
    cmake_parse_arguments("KAGEN"
        ""
        ""
        "FILES;CORES"
        ${ARGN})

    katestrophe_add_test_executable(${KAGEN_TARGET_NAME} FILES ${KAGEN_FILES})
    target_link_libraries(${KAGEN_TARGET_NAME} PRIVATE
        mpi-gtest-main
        KaGen::KaGen)
    target_include_directories(${KAGEN_TARGET_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    katestrophe_add_mpi_test(${KAGEN_TARGET_NAME} CORES ${KAGEN_CORES})
endfunction ()

###
### Add unit tests
###

# Path Graph Generator
kagen_add_test(test_path_generator
    FILES path/general_path_generator_test.cpp
    CORES 1 2 3 4 8 16)

# Image Mesh Generator
kagen_add_test(test_kargb 
    FILES image/kargb_test.cpp 
    CORES 1)

# File Graph Generator
kagen_add_test(test_generic_file_generator
    FILES file/generic_file_generator_test.cpp
    CORES 1 2 3 4 8 16)

# Geometric Generator
kagen_add_test(test_rgg2d
        FILES geometric/rgg2d_test.cpp
        CORES 1 2 3 4 5 6)

kagen_add_test(test_rgg3d
        FILES geometric/rgg3d_test.cpp
        CORES 1 2 3 4 5 6)

kagen_add_test(test_hyperbolic
        FILES geometric/hyperbolic_test.cpp
        CORES 1 2 3 4 5 6)

