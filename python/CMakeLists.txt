find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(_kagen src/_kagen.cpp)
target_link_libraries(_kagen PRIVATE KaGen)
target_include_directories(_kagen PRIVATE ${PROJECT_SOURCE_DIR})

install(TARGETS _kagen DESTINATION kagen)
