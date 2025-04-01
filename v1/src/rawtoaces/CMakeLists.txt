cmake_minimum_required(VERSION 3.10)
include_directories( "${CMAKE_CURRENT_SOURCE_DIR}" )

add_executable( rawtoaces_v1
    main.cpp
)

set_property(TARGET rawtoaces_v1 PROPERTY CXX_STANDARD 17)

target_link_libraries ( rawtoaces_v1
    PUBLIC
        ${RAWTOACESLIB}
)

if ( LIBRAW_CONFIG_FOUND )
    target_link_libraries ( rawtoaces_v1 PUBLIC libraw::raw )
else ()
    target_include_directories(rawtoaces_v1 PUBLIC ${libraw_INCLUDE_DIR}/.. )    
    target_link_directories(rawtoaces_v1 PUBLIC ${libraw_LIBRARY_DIRS} )
    target_link_libraries(rawtoaces_v1 PUBLIC ${libraw_LIBRARIES} ${libraw_LDFLAGS_OTHER} )
endif ()

install( TARGETS rawtoaces_v1 DESTINATION bin )
