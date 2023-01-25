include(FindPkgConfig)
pkg_check_modules(ImageMagick IMPORTED_TARGET ImageMagick Magick++ MagickCore)

if (ImageMagick_FOUND)
    add_executable(img2kargb img2kargb.cpp)
    target_link_libraries(img2kargb PRIVATE PkgConfig::ImageMagick)
    install(TARGETS img2kargb DESTINATION bin)

    add_executable(part2img part2img.cpp)
    target_link_libraries(part2img PRIVATE PkgConfig::ImageMagick)
    install(TARGETS part2img DESTINATION bin)
else ()
    message(STATUS "img2kargb not available: could not find ImageMagick on your system")
    message(STATUS "part2img not available: could not find ImageMagick on your system")
endif ()

add_executable(upsb2kargb upsb2kargb.cpp)
install(TARGETS upsb2kargb DESTINATION bin)
