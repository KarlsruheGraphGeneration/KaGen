################################################################################
# app/CMakeLists.txt
#
# CMake file for generator app
#
# Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
#
# All rights reserved. Published under the BSD-2 license in the LICENSE file.
################################################################################

add_executable(KaGenApp KaGen.cpp)
set_property(TARGET KaGenApp PROPERTY OUTPUT_NAME KaGen)
target_link_libraries(KaGenApp PRIVATE KaGen::KaGen)
target_compile_features(KaGenApp PRIVATE cxx_std_17)

if (KAGEN_BUILD_TOOLS)
    add_subdirectory(tools)
endif ()
