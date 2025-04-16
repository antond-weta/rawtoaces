cmake_minimum_required(VERSION 3.10)
include_directories( "${CMAKE_CURRENT_SOURCE_DIR}" )

add_library ( rawtoaces_util ${DO_SHARED}
    data_adapters/spectral_loader.cpp
    data_adapters/spectral_loader.hpp
    transform_cache.cpp
    transform_cache.hpp
    rawtoaces_util.cpp
    rawtoaces_util.hpp
    cache_base.hpp
)

set_property(TARGET rawtoaces_util PROPERTY CXX_STANDARD 20)

if ( OIIO_FOUND )
    target_include_directories ( rawtoaces_util PRIVATE ${OIIO_INCLUDE_DIRS} )
    target_link_directories    ( rawtoaces_util PUBLIC  ${OIIO_LIBRARY_DIRS} )
    target_link_libraries      ( rawtoaces_util
        PUBLIC
            ${OIIO_LIBRARIES}
            ${OIIO_LDFLAGS_OTHER}
    )
endif()
 
target_link_libraries ( rawtoaces_util
    PRIVATE
        nlohmann_json::nlohmann_json
    PUBLIC
        rawtoaces_core
        OpenImageIO::OpenImageIO
)
    
set_target_properties( rawtoaces_util PROPERTIES
    SOVERSION ${RTA_MAJOR_VERSION}.${RTA_MINOR_VERSION}.${RTA_PATCH_VERSION}
    VERSION ${RTA_VERSION} )

install(
    FILES
        ${PROJECT_SOURCE_DIR}/include/rawtoaces/rawtoaces_util.h
        ${PROJECT_SOURCE_DIR}/include/rawtoaces/cache_base.h	 	
    DESTINATION include/rawtoaces
)

install( TARGETS rawtoaces_util DESTINATION lib )
