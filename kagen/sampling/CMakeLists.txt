################################################################################
# sampling/CMakeLists.txt
#
# CMake file for sampling library
#
# Copyright (C) 2017 Lorenz Hübschle-Schneider <lorenz@4z2.de>
#
# All rights reserved. Published under the BSD-2 license in the LICENSE file.
################################################################################

# glob general sources
file(GLOB SAMPLING_SRCS
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  *.[ch]pp
  rng/stl.hpp
  rng/dSFMT.hpp
  rng/dSFMT.cpp
  )

# if we have Intel MKL, add rng/mkl.hpp
if(MKL_FOUND)
  set(SAMPLING_SRCS ${SAMPLING_SRCS} rng/mkl.hpp)
endif()

add_library(sampling STATIC ${SAMPLING_SRCS})
set_target_properties(sampling PROPERTIES LINKER_LANGUAGE CXX)
target_compile_definitions(sampling PUBLIC ${SAMPLING_DEFINITIONS})
target_include_directories(sampling PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(sampling SYSTEM PUBLIC ${PROJECT_SOURCE_DIR}/extlib/spooky)
target_include_directories(sampling SYSTEM PUBLIC ${SAMPLING_INCLUDE_DIRS})
target_link_libraries(sampling ${SAMPLING_LINK_LIBRARIES})

################################################################################
