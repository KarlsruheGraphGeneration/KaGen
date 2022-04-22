################################################################################
# app/CMakeLists.txt
#
# CMake file for generator app
#
# Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
#
# All rights reserved. Published under the BSD-2 license in the LICENSE file.
################################################################################

add_executable(kagen kagen.cpp)
target_link_libraries(kagen PRIVATE kagen_obj)
target_compile_features(kagen INTERFACE cxx_std_17)

################################################################################
