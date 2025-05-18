cmake_minimum_required(VERSION 3.10)
include_directories( "${CMAKE_CURRENT_SOURCE_DIR}" )

add_executable( rawtoaces
    main.cpp
)

set_property(TARGET rawtoaces PROPERTY CXX_STANDARD 20)

target_include_directories( rawtoaces
    PUBLIC
        ${OpenImageIO_INCLUDE_DIR}
)

target_link_libraries ( rawtoaces
    PUBLIC
#        OpenImageIO::OpenImageIO
        rawtoaces::util
)


set_target_properties(rawtoaces PROPERTIES XCODE_SCHEME_ARGUMENTS "${CMAKE_XCODE_SCHEME_ARGUMENTS}" )

install( TARGETS rawtoaces DESTINATION bin )
