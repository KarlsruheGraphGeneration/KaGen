include(FindPkgConfig)
pkg_check_modules(ImageMagick IMPORTED_TARGET ImageMagick Magick++ MagickCore)

if (ImageMagick_FOUND AND ImageMagick_ImageMagick_VERSION VERSION_GREATER_EQUAL "7.0.0")
    add_executable(img2kargb img2kargb.cpp)
    target_link_libraries(img2kargb PRIVATE PkgConfig::ImageMagick)
    target_include_directories(img2kargb PRIVATE ${PROJECT_SOURCE_DIR})
    install(TARGETS img2kargb DESTINATION bin)

    add_executable(part2img part2img.cpp)
    target_link_libraries(part2img PRIVATE PkgConfig::ImageMagick)
    target_include_directories(part2img PRIVATE ${PROJECT_SOURCE_DIR})
    install(TARGETS part2img 
        COMPONENT tools
        EXCLUDE_FROM_ALL
        DESTINATION bin)
else ()
    if (NOT ImageMagick_FOUND) 
        set(REASON "could not find ImageMagick on your system")
    else ()
        set(REASON "ImageMagick found on your system is too old (version ${ImageMagick_ImageMagick_VERSION}, require >= 7.0.0)")
    endif ()
    message(STATUS "img2kargb not available: ${REASON}")
    message(STATUS "part2img not available: ${REASON}")
endif ()

add_executable(upsb2kargb upsb2kargb.cpp)
target_include_directories(upsb2kargb PRIVATE ${PROJECT_SOURCE_DIR})

add_executable(chkgraph chkgraph.cpp)
target_link_libraries(chkgraph PRIVATE kagen)
target_include_directories(chkgraph PRIVATE ${PROJECT_SOURCE_DIR})

add_executable(graphstats graphstats.cpp)
target_link_libraries(graphstats PRIVATE kagen)
target_include_directories(graphstats PRIVATE ${PROJECT_SOURCE_DIR})

add_executable(pangraph pangraph.cpp)
target_link_libraries(pangraph PRIVATE kagen)
target_include_directories(pangraph PRIVATE ${PROJECT_SOURCE_DIR})

install(TARGETS upsb2kargb graphstats chkgraph pangraph
    COMPONENT tools
    EXCLUDE_FROM_ALL
    DESTINATION bin)
