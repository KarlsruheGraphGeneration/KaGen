add_executable(example_csr csr_example.cpp)
target_link_libraries(example_csr PRIVATE KaGen::KaGen MPI::MPI_CXX)
