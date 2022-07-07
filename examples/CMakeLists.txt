add_executable(example_csr csr_example.cpp)
target_link_libraries(example_csr PRIVATE KaGen::KaGen)

add_executable(example_rgg2d rgg2d_example.cpp)
target_link_libraries(example_rgg2d PRIVATE KaGen::KaGen)
