cmake_minimum_required(VERSION 3.10)
include_directories( "${CMAKE_CURRENT_SOURCE_DIR}" )


add_compile_options("$<$<C_COMPILER_ID:MSVC>:/utf-8>")
add_compile_options("$<$<CXX_COMPILER_ID:MSVC>:/utf-8>")

enable_testing()

################################################################################

add_executable( test_core_math
    test_data.hpp
    core_math.cpp
)

target_link_libraries(
    test_core_math
    PUBLIC
        rawtoaces::core
        OpenImageIO::OpenImageIO
)

set_property(TARGET test_core_math PROPERTY CXX_STANDARD 20)
add_test ( NAME test_core_math COMMAND test_core_math )

################################################################################

add_executable( test_core_spectral
    test_data.hpp
    core_spectral.cpp
)

target_link_libraries(
    test_core_spectral
    PUBLIC
        rawtoaces::core
        OpenImageIO::OpenImageIO
)

set_property(TARGET test_core_spectral PROPERTY CXX_STANDARD 20)
add_test ( NAME test_core_spectral COMMAND test_core_spectral )

################################################################################

add_executable( test_core_metadata
    test_data.hpp
    core_metadata.cpp
)

target_link_libraries(
    test_core_metadata
    PUBLIC
        rawtoaces::core
        OpenImageIO::OpenImageIO
)

set_property(TARGET test_core_metadata PROPERTY CXX_STANDARD 20)
add_test ( NAME test_core_metadata COMMAND test_core_metadata )

################################################################################

add_executable( test_util_spectral
    test_data.hpp
    util_spectral.cpp
)

target_link_libraries(
    test_util_spectral
    PUBLIC
        rawtoaces::core
        rawtoaces::util
        OpenImageIO::OpenImageIO
)

set_property(TARGET test_util_spectral PROPERTY CXX_STANDARD 20)
add_test ( NAME test_util_spectral COMMAND test_util_spectral )

################################################################################



if (${RTA_ENABLE_EXIFTOOL})

add_executable( test_util_exiftool
    util_exiftool.cpp
)

if (CMAKE_GENERATOR STREQUAL "Xcode")
    target_compile_definitions(test_util_exiftool PUBLIC RTA_XCODE)
endif (CMAKE_GENERATOR STREQUAL "Xcode")

target_link_libraries(
    test_util_exiftool
    PUBLIC
        rawtoaces::util
        OpenImageIO::OpenImageIO
)

set_property(TARGET test_util_exiftool PROPERTY CXX_STANDARD 20)
add_test ( NAME test_util_exiftool COMMAND test_util_exiftool )

endif (${RTA_ENABLE_EXIFTOOL})
