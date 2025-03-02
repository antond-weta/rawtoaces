cmake_minimum_required(VERSION 3.10)
include_directories( "${CMAKE_CURRENT_SOURCE_DIR}" )

add_library( rawtoaces_core ${DO_SHARED}
    spectral_solver.cpp
    matrix_solver.cpp

    # Make the headers visible in IDEs. This should not affect the builds.
    ../../include/rawtoaces/define.h
    ../../include/rawtoaces/mathOps.h
    ../../include/rawtoaces/rta.h
    ../../include/rawtoaces/metadata.h

    ../../include/rawtoaces/spectrum.hpp
    ../../include/rawtoaces/spectral_data.hpp
    ../../include/rawtoaces/spectral_solver.hpp
    ../../include/rawtoaces/matrix_solver.hpp
)

set_property(TARGET rawtoaces_core PROPERTY CXX_STANDARD 20)

# Eigen is also used in rawtoaces_util and rawtoaces_util2, so setting it public here
target_link_libraries(
    rawtoaces_core
    PUBLIC
        Eigen3::Eigen
    PRIVATE
        nlohmann_json::nlohmann_json
)

# Ceres is also used in rawtoaces_util and rawtoaces_util2, so setting it public here
if ( ${Ceres_VERSION_MAJOR} GREATER 1 )
    target_link_libraries( rawtoaces_core PUBLIC Ceres::ceres )
else ()
    target_include_directories(rawtoaces_core PUBLIC ${CERES_INCLUDE_DIRS})
    target_link_libraries(rawtoaces_core PUBLIC ${CERES_LIBRARIES})
endif ()


set_target_properties( rawtoaces_core PROPERTIES
    SOVERSION ${RAWTOACES_MAJOR_VERSION}.${RAWTOACES_MINOR_VERSION}.${RAWTOACES_PATCH_VERSION}
    VERSION ${RAWTOACES_VERSION}
)

install(FILES
    ${PROJECT_SOURCE_DIR}/include/rawtoaces/define.h
    ${PROJECT_SOURCE_DIR}/include/rawtoaces/mathOps.h
    ${PROJECT_SOURCE_DIR}/include/rawtoaces/rta.h
    ${PROJECT_SOURCE_DIR}/include/rawtoaces/metadata.h
    ${PROJECT_SOURCE_DIR}/include/rawtoaces/spectrum
    ${PROJECT_SOURCE_DIR}/include/rawtoaces/spectral_data.h
    ${PROJECT_SOURCE_DIR}/include/rawtoaces/spectral_solver.h
    ${PROJECT_SOURCE_DIR}/include/rawtoaces/matrix_solver.h
    DESTINATION include/rawtoaces
)

install( TARGETS rawtoaces_core DESTINATION lib )
