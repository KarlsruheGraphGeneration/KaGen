find_package(ImageMagick QUIET COMPONENTS Magick++ MagickCore)

if (ImageMagick_FOUND AND ImageMagick_Magick++_FOUND)
    # @todo is there a way to deduce these values from the library name found by FindImageMagick.cmake?
    set(KAGEN_IM_DEPTH 16)
    set(KAGEN_IM_HDRI 1)

    add_executable(img2kargb img2kargb.cpp)
    target_compile_options(img2kargb PRIVATE -DMAGICKCORE_QUANTUM_DEPTH=${KAGEN_IM_DEPTH} -DMAGICKCORE_HDRI_ENABLE=${KAGEN_IM_HDRI})
    target_link_libraries(img2kargb PRIVATE ${ImageMagick_LIBRARIES} ${ImageMagick_Magick++_LIBRARIES} ${ImageMagick_MagickCore_LIBRARIES})
    target_include_directories(img2kargb SYSTEM PRIVATE ${ImageMagick_INCLUDE_DIRS} ${ImageMagick_Magick++_INCLUDE_DIRS} ${ImageMagick_MagickCore_INCLUDE_DIRS})
    install(TARGETS img2kargb DESTINATION bin)

    add_executable(part2img part2img.cpp)
    target_compile_options(part2img PRIVATE -DMAGICKCORE_QUANTUM_DEPTH=${KAGEN_IM_DEPTH} -DMAGICKCORE_HDRI_ENABLE=${KAGEN_IM_HDRI})
    target_link_libraries(part2img PRIVATE ${ImageMagick_LIBRARIES} ${ImageMagick_Magick++_LIBRARIES} ${ImageMagick_MagickCore_LIBRARIES})
    target_include_directories(part2img SYSTEM PRIVATE ${ImageMagick_INCLUDE_DIRS} ${ImageMagick_Magick++_INCLUDE_DIRS} ${ImageMagick_MagickCore_INCLUDE_DIRS})
    install(TARGETS part2img DESTINATION bin)
else ()
    message(STATUS "img2kargb not available: could not find ImageMagick on your system")
    message(STATIS "part2img not available: could not find ImageMagick on your system")
endif ()
