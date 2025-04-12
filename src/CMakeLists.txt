cmake_minimum_required(VERSION 3.10)

add_compile_options("$<$<C_COMPILER_ID:MSVC>:/utf-8>")
add_compile_options("$<$<CXX_COMPILER_ID:MSVC>:/utf-8>")

add_subdirectory("core")

set(INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include/
    CACHE PATH "Location of header files" )

install( EXPORT rawtoacesTargets
    FILE rawtoacesTargets.cmake
    DESTINATION lib/cmake/rawtoaces
)

include(CMakePackageConfigHelpers)

configure_package_config_file(
    Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/rawtoacesConfig.cmake
    INSTALL_DESTINATION lib/cmake/rawtoaces
    PATH_VARS INCLUDE_INSTALL_DIR
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/rawtoacesConfigVersion.cmake
    VERSION ${RTA_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/rawtoacesConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/rawtoacesConfigVersion.cmake
    DESTINATION lib/cmake/rawtoaces
)
