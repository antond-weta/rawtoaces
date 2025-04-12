cmake_minimum_required(VERSION 3.10)
include_directories( "${CMAKE_CURRENT_SOURCE_DIR}" )

add_library( rawtoaces_core ${DO_SHARED}
    constants.hpp
    spectrum.hpp
    metadata.hpp
    spectral_data.hpp

    metadata_solver.hpp
    metadata_solver.cpp

    spectral_solver.hpp
    spectral_solver.cpp

    matrix_solver.hpp
    matrix_solver.cpp
)

set_property(TARGET rawtoaces_core PROPERTY CXX_STANDARD 20)

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

set_target_properties( rawtoaces_core PROPERTIES
    SOVERSION ${RTA_MAJOR_VERSION}.${RTA_MINOR_VERSION}.${RTA_PATCH_VERSION}
    VERSION ${RTA_VERSION}
)

install(FILES
    constants.hpp
    spectrum.hpp
    metadata.hpp
    metadata_solver.hpp
    spectral_solver.hpp
    matrix_solver.hpp
    DESTINATION include/rawtoaces
)

install( TARGETS rawtoaces_core DESTINATION lib )
