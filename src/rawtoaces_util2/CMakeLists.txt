# Copyright Contributors to the rawtoaces project.
# SPDX-License-Identifier: Apache-2.0
# https://github.com/AcademySoftwareFoundation/rawtoaces

include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")

add_library (${RawToAcesUtil2_LIB} ${DO_SHARED}
    transform_cache.cpp
    transform_cache.h
    image_converter.cpp
    ${RawToAces_INCLUDE_DIR}/rawtoaces/image_converter.h
    ${RawToAces_INCLUDE_DIR}/rawtoaces/cache_base.h
)

set_property (TARGET ${RawToAcesUtil2_LIB} PROPERTY CXX_STANDARD 17)

if (OIIO_FOUND)
    target_include_directories (${RawToAcesUtil2_LIB} PRIVATE ${OIIO_INCLUDE_DIRS})
    target_link_directories (${RawToAcesUtil2_LIB} PUBLIC  ${OIIO_LIBRARY_DIRS})
    target_link_libraries (${RawToAcesUtil2_LIB}
        PUBLIC
            ${OIIO_LIBRARIES}
            ${OIIO_LDFLAGS_OTHER}
    )
endif ()
 
target_link_libraries (${RawToAcesUtil2_LIB}
    PRIVATE
        ${RawToAcesIdt_LIB}
    PUBLIC
        OpenImageIO::OpenImageIO
)
    
set_target_properties(
    ${RawToAcesUtil2_LIB} PROPERTIES
    SOVERSION ${RawToAces_MAJOR_VERSION}
    VERSION ${RawToAces_VERSION}
)

install(
    FILES
         ${RawToAces_INCLUDE_DIR}/rawtoaces/image_converter.h
         ${RawToAces_INCLUDE_DIR}/rawtoaces/cache_base.h	 	
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
)

install (TARGETS ${RawToAcesUtil2_LIB} DESTINATION ${CMAKE_INSTALL_LIBDIR})
