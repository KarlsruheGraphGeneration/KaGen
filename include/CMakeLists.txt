################################################################################
# include/CMakeLists.txt
#
# CMake file for generator library
#
# Copyright (C) 2017 Sebastian Lamm <lamm@kit.edu>
#
# All rights reserved. Published under the BSD-2 license in the LICENSE file.
################################################################################

# glob general sources
file(GLOB GENERATOR_SRCS
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*.[ch]pp
  )

add_library(generators STATIC ${GENERATOR_SRCS})
set_target_properties(generators PROPERTIES LINKER_LANGUAGE CXX)
target_compile_definitions(generators PUBLIC ${KAGEN_DEFINITIONS})
target_include_directories(generators PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(generators PUBLIC ${PROJECT_SOURCE_DIR}/extlib/spooky)
target_include_directories(generators SYSTEM PUBLIC ${KAGEN_INCLUDE_DIRS})
target_link_libraries(generators ${KAGEN_LINK_LIBRARIES})

################################################################################
