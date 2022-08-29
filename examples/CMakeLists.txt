add_executable(example_csr csr_example.cpp)
target_link_libraries(example_csr PRIVATE KaGen::KaGen)

add_executable(example_c_csr csr_c_example.c)
target_link_libraries(example_c_csr PRIVATE KaGen::cKaGen)

add_executable(example_rgg2d rgg2d_example.cpp)
target_link_libraries(example_rgg2d PRIVATE KaGen::KaGen)

add_executable(example_c_rgg2d rgg2d_c_example.c)
target_link_libraries(example_c_rgg2d PRIVATE KaGen::cKaGen)
