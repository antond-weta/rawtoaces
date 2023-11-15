# Until we get some of these modules into the upstream packages, put them here
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_INSTALL_PREFIX}/share/CMake")


find_package ( Eigen3        CONFIG REQUIRED )
# find_package ( Imath         CONFIG REQUIRED )
find_package ( Ceres                REQUIRED )
find_package ( Boost                REQUIRED
    COMPONENTS
        system
        filesystem
        unit_test_framework
)

if ( USE_OPENEXR )
   # find_package ( OpenEXR CONFIG QUIET )
    
    if ( OpenEXR_FOUND )
        message ( STATUS "OpenEXR config found ${openexr_VERSION}" )
        set ( OPENEXR_CONFIG_FOUND TRUE )
    else ()
        message ( WARNING "OpenEXR config not found" )
        find_package ( OpenEXR MODULE REQUIRED )
    endif ()
else ()
    find_package ( AcesContainer CONFIG REQUIRED )
endif ()

find_package (libraw CONFIG QUIET )

if (libraw_FOUND )
    message("STATUS LibRaw config found")
    set ( LIBRAW_CONFIG_FOUND TRUE )
else ()
    message("WARNING LibRaw config not found, trying to find a module.")
    find_package(libraw MODULE REQUIRED)
endif ()
