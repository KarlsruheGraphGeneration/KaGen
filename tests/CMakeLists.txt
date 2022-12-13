function (add_kagen_test target)
    add_executable(${target} ${ARGN})
    target_link_libraries(${target} PRIVATE GTest::gtest_main kagen_obj)
    gtest_discover_tests(${target})
endfunction ()

add_kagen_test(test_image_mesh image_mesh_test.cpp)

