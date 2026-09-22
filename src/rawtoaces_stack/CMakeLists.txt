cmake_minimum_required(VERSION 3.12)

# Include coverage support if enabled
if( ENABLE_COVERAGE )
    include( ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/CodeCoverage.cmake )
endif()

add_executable( rawtoaces_stack
    main.cpp
)

target_link_libraries ( rawtoaces_stack
    PUBLIC
        ${RAWTOACES_UTIL_LIB}
)

# Enable coverage for this executable if coverage is enabled
if( ENABLE_COVERAGE AND COVERAGE_SUPPORTED )
    setup_coverage_flags(rawtoaces_stack)
endif()

install( TARGETS rawtoaces_stack DESTINATION ${INSTALL_BIN_DIR} )

# Simplifies debugging with Xcode, doesn't affect anything else.
set_target_properties(rawtoaces_stack PROPERTIES XCODE_SCHEME_ARGUMENTS "${CMAKE_XCODE_SCHEME_ARGUMENTS}" )