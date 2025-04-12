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
        rawtoaces_core
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
        rawtoaces_core
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
        rawtoaces_core
        OpenImageIO::OpenImageIO
)

set_property(TARGET test_core_metadata PROPERTY CXX_STANDARD 20)
add_test ( NAME test_core_metadata COMMAND test_core_metadata )

