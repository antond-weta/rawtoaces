cmake_minimum_required(VERSION 3.10)
include_directories( "${CMAKE_CURRENT_SOURCE_DIR}" )


add_compile_options("$<$<C_COMPILER_ID:MSVC>:/utf-8>")
add_compile_options("$<$<CXX_COMPILER_ID:MSVC>:/utf-8>")

enable_testing()

################################################################################

add_executable( test_core
    core_test.cpp
)

target_link_libraries(
    test_core
    PUBLIC
        rawtoaces::core
        OpenImageIO::OpenImageIO
)

set_property(TARGET test_core PROPERTY CXX_STANDARD 20)
add_test ( NAME test_core COMMAND test_core )
