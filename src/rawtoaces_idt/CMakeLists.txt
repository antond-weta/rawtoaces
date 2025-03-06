# Copyright Contributors to the rawtoaces project.
# SPDX-License-Identifier: Apache-2.0
# https://github.com/AcademySoftwareFoundation/rawtoaces

include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

add_library (
    ${RawToAcesIdt_LIB} ${DO_SHARED}
    rta.cpp
    # make the headers visible in IDEs. This should not affect the builds.
    ${RawToAces_INCLUDE_DIR}/rawtoaces/define.h
    ${RawToAces_INCLUDE_DIR}/rawtoaces/mathops.h
    ${RawToAces_INCLUDE_DIR}/rawtoaces/rta.h
    ${RawToAces_INCLUDE_DIR}/rawtoaces/metadata.h
)
set_property (TARGET ${RawToAcesIdt_LIB} PROPERTY CXX_STANDARD 17)

# eigen is also used in rawtoaces_util and rawtoaces_util2, so setting it public here
target_link_libraries (
    ${RawToAcesIdt_LIB}
    PUBLIC
        Eigen3::Eigen
    PRIVATE
	    nlohmann_json::nlohmann_json
)

# ceres is also used in rawtoaces_util and rawtoaces_util2, so setting it public here
if (${Ceres_VERSION_MAJOR} GREATER 1 )
    target_link_libraries (${RawToAcesIdt_LIB} PUBLIC Ceres::ceres)
else ()
    target_include_directories (${RawToAcesIdt_LIB} PUBLIC ${CERES_INCLUDE_DIRS})
    target_link_libraries (${RawToAcesIdt_LIB} PUBLIC ${CERES_LIBRARIES})
endif ()

set_target_properties(
    ${RawToAcesIdt_LIB} PROPERTIES
    SOVERSION ${RawToAces_MAJOR_VERSION}
    VERSION ${RawToAces_VERSION}
)

install (FILES
    ${RawToAces_INCLUDE_DIR}/rawtoaces/define.h
    ${RawToAces_INCLUDE_DIR}/rawtoaces/mathops.h
    ${RawToAces_INCLUDE_DIR}/rawtoaces/rta.h
    ${RawToAces_INCLUDE_DIR}/rawtoaces/metadata.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/rawtoaces
)

install (TARGETS ${RawToAcesIdt_LIB} DESTINATION ${CMAKE_INSTALL_LIBDIR})
