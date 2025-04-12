cmake_minimum_required(VERSION 3.10)
#include_directories( "${CMAKE_CURRENT_SOURCE_DIR}" )

set( CORE_PUBLIC_HEADER
    constants.hpp
    spectrum.hpp
    metadata.hpp
    spectral_data.hpp
    metadata_solver.hpp
    spectral_solver.hpp
    matrix_solver.hpp
)

add_library( rawtoaces_core ${DO_SHARED}
    ${CORE_PUBLIC_HEADER}
    metadata_solver.cpp
    spectral_solver.cpp
    matrix_solver.cpp
)

add_library( rawtoaces::core ALIAS rawtoaces_core )

set_target_properties( rawtoaces_core PROPERTIES
    CXX_STANDARD 20
    OUTPUT_NAME "rawtoaces_core"
    EXPORT_NAME "rawtoaces::core"
    PUBLIC_HEADER "${CORE_PUBLIC_HEADER}"
    SOVERSION ${RTA_MAJOR_VERSION}.${RTA_MINOR_VERSION}.${RTA_PATCH_VERSION}
    VERSION ${RTA_VERSION}
)

target_include_directories( rawtoaces_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>
    $<INSTALL_INTERFACE:./include>
)

if( RTA_ENABLE_EIGEN )
    target_link_libraries( rawtoaces_core PUBLIC Eigen3::Eigen )
    target_compile_definitions( rawtoaces_core PRIVATE ENABLE_EIGEN=1 )
endif( RTA_ENABLE_EIGEN )

if( RTA_ENABLE_CERES )
    if ( ${Ceres_VERSION_MAJOR} GREATER 1 )
        target_link_libraries( rawtoaces_core PUBLIC Ceres::ceres )
    else ()
        target_include_directories( rawtoaces_core PUBLIC ${CERES_INCLUDE_DIRS})
        target_link_libraries( rawtoaces_core PUBLIC ${CERES_LIBRARIES})
    endif ()

    target_compile_definitions( rawtoaces_core PRIVATE ENABLE_CERES=1 )
endif( RTA_ENABLE_CERES )

install( TARGETS rawtoaces_core
    EXPORT rawtoacesTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    #ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    #RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    #LIBRARY DESTINATION lib
    PUBLIC_HEADER DESTINATION include/rawtoaces/core
)
