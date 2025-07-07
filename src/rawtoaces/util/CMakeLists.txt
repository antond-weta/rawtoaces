cmake_minimum_required( VERSION 3.10 )
#include_directories( "${CMAKE_CURRENT_SOURCE_DIR}" )

set( UTIL_PUBLIC_HEADER
    rawtoaces_util.hpp
)

add_library ( rawtoaces_util ${DO_SHARED}
    ${UTIL_PUBLIC_HEADER}

    data_adapters/spectral_loader.cpp
    data_adapters/spectral_loader.hpp
    transform_cache.cpp
    transform_cache.hpp

    cache_base.hpp

    $<$<BOOL:${RTA_ENABLE_EXIFTOOL}>:
        data_adapters/exiftool.cpp
        data_adapters/exiftool.hpp
    >

    rawtoaces_util.cpp
)

add_library ( rawtoaces::util ALIAS rawtoaces_util )

set_target_properties( rawtoaces_util PROPERTIES
    CXX_STANDARD 20
    OUTPUT_NAME "rawtoaces_util"
    EXPORT_NAME "rawtoaces::util"
    PUBLIC_HEADER "${UTIL_PUBLIC_HEADER}"
    SOVERSION ${RTA_MAJOR_VERSION}.${RTA_MINOR_VERSION}.${RTA_PATCH_VERSION}
    VERSION ${RTA_VERSION} 
)

target_include_directories( rawtoaces_util PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>
    $<INSTALL_INTERFACE:./include>
)

#target_include_directories ( rawtoaces_util PRIVATE ${OIIO_INCLUDE_DIRS} )
#target_link_directories    ( rawtoaces_util PUBLIC  ${OIIO_LIBRARY_DIRS} )
#target_link_libraries      ( rawtoaces_util
#    PUBLIC
#        ${OIIO_LIBRARIES}
#        ${OIIO_LDFLAGS_OTHER}
#)
 
target_link_libraries ( rawtoaces_util
    PRIVATE
        nlohmann_json::nlohmann_json
    PUBLIC
        rawtoaces::core
        OpenImageIO::OpenImageIO
)

target_compile_definitions( rawtoaces_util PRIVATE RTA_VERSION="${RTA_VERSION}" )

if ( ${RTA_ENABLE_EXIFTOOL} )
    target_compile_definitions(rawtoaces_util PRIVATE RTA_ENABLE_EXIFTOOL=1 )
endif()

install( TARGETS rawtoaces_util
    EXPORT rawtoacesTargets
    LIBRARY DESTINATION lib
    PUBLIC_HEADER DESTINATION include/rawtoaces/util
)
