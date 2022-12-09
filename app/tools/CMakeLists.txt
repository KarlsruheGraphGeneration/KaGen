find_package(ImageMagick QUIET COMPONENTS Magick++ MagickCore)

if (ImageMagick_FOUND AND ImageMagick_Magick++_FOUND)
    add_executable(img2kargb img2kargb.cpp)
    target_compile_options(img2kargb PRIVATE -DMAGICKCORE_QUANTUM_DEPTH=16 -DMAGICKCORE_HDRI_ENABLE=1)
    target_link_libraries(img2kargb PRIVATE ${ImageMagick_LIBRARIES} ${ImageMagick_Magick++_LIBRARIES} ${ImageMagick_MagickCore_LIBRARIES})
    target_include_directories(img2kargb SYSTEM PRIVATE ${ImageMagick_INCLUDE_DIRS} ${ImageMagick_Magick++_INCLUDE_DIRS} ${ImageMagick_MagickCore_INCLUDE_DIRS})

    install(TARGETS img2kargb DESTINATION bin)
else ()
    message(STATUS "img2kargb not available: could not find ImageMagick on your system")
endif ()
