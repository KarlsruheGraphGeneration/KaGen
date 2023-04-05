add_executable(example_cpp_csr csr_cpp_example.cpp)
target_link_libraries(example_cpp_csr PRIVATE KaGen::KaGen)

add_executable(example_c_csr csr_c_example.c)
target_link_libraries(example_c_csr PRIVATE KaGen::KaGen)

add_executable(example_cpp_rgg2d rgg2d_cpp_example.cpp)
target_link_libraries(example_cpp_rgg2d PRIVATE KaGen::KaGen)

add_executable(example_c_rgg2d rgg2d_c_example.c)
target_link_libraries(example_c_rgg2d PRIVATE KaGen::KaGen)

add_executable(example_cpp_options_string options_string_cpp_example.cpp)
target_link_libraries(example_cpp_options_string PRIVATE KaGen::KaGen)

add_executable(example_c_options_string options_string_c_example.c)
target_link_libraries(example_c_options_string PRIVATE KaGen::KaGen)

