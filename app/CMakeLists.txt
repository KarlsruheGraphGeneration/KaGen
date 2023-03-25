################################################################################
# app/CMakeLists.txt
#
# CMake file for generator app
#
# Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
#
# All rights reserved. Published under the BSD-2 license in the LICENSE file.
################################################################################

add_executable(KaGen KaGen.cpp)
target_link_libraries(KaGen PRIVATE kagen)
target_compile_features(KaGen PRIVATE cxx_std_17)

install(TARGETS KaGen DESTINATION bin)

add_subdirectory(tools)
