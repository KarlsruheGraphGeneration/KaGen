add_library(kagen kagen.cpp kagen.h)
target_link_libraries(kagen PUBLIC kagen_obj)
target_include_directories(kagen PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
