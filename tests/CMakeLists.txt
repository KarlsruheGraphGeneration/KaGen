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
        kagen_obj
        KaGen::KaGen)
    target_include_directories(${KAGEN_TARGET_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    katestrophe_add_mpi_test(${KAGEN_TARGET_NAME} CORES ${KAGEN_CORES})
endfunction ()

###
### Add unit tests
###

# Image Mesh Generator
kagen_add_test(test_kargb 
    FILES image/kargb_test.cpp 
    CORES 1)

# Static Graph Generator
kagen_add_test(test_metis_generator
    FILES static/metis_test.cpp
    CORES 1 2 3 4 8 16)

#kagen_add_test(test_binary_parhip_generator
#FILES static/binary_parhip_test.cpp
#CORES 1 2 3 4 8 16)

