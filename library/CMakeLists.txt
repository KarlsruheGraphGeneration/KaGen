add_library(kagen_library kagen_library.cpp kagen_library.h)
target_link_libraries(kagen_library PUBLIC kagen_obj)
target_include_directories(kagen_library PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
