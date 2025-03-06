# Copyright Contributors to the rawtoaces project.
# SPDX-License-Identifier: Apache-2.0
# https://github.com/AcademySoftwareFoundation/rawtoaces

include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")

add_library (
    ${RawToAcesUtil_LIB} ${DO_SHARED}
    acesrender.cpp
    # make the headers visible in IDEs. This should not affect the builds.
    ${RawToAces_INCLUDE_DIR}/rawtoaces/acesrender.h
)

set_property(TARGET ${RawToAcesUtil_LIB} PROPERTY CXX_STANDARD 17)

target_include_directories (${RawToAcesUtil_LIB} PRIVATE ${AcesContainer_INCLUDE_DIRS})
target_link_directories (${RawToAcesUtil_LIB} PUBLIC ${AcesContainer_LIBRARY_DIRS})
target_link_libraries (${RawToAcesUtil_LIB}
    PUBLIC
        ${AcesContainer_LIBRARIES}
        ${AcesContainer_LDFLAGS_OTHER}
)
 
target_link_libraries (${RawToAcesUtil_LIB}
    PUBLIC
        ${RawToAcesIdt_LIB}
    PRIVATE
        Imath::Imath
        Imath::ImathConfig
)
 
if (LIBRAW_CONFIG_FOUND)
    target_link_libraries (${RawToAcesUtil_LIB} PUBLIC libraw::raw)
else ()
    target_include_directories (${RawToAcesUtil_LIB} PUBLIC ${libraw_INCLUDE_DIR}/..)
    target_link_directories (${RawToAcesUtil_LIB} PUBLIC ${libraw_LIBRARY_DIRS})
    target_link_libraries (${RawToAcesUtil_LIB} PUBLIC ${libraw_LIBRARIES} ${libraw_LDFLAGS_OTHER})
endif ()

set_target_properties(
    ${RawToAcesUtil_LIB} PROPERTIES
    SOVERSION ${RawToAces_MAJOR_VERSION}
    VERSION ${RawToAces_VERSION}
)
       
install(FILES
    ${RawToAces_INCLUDE_DIR}/rawtoaces/acesrender.h	 	
 	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
)

install (TARGETS ${RAWTOACESLIB} DESTINATION ${CMAKE_INSTALL_LIBDIR})
