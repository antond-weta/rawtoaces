# Copyright Contributors to the rawtoaces project.
# SPDX-License-Identifier: Apache-2.0
# https://github.com/AcademySoftwareFoundation/rawtoaces

cmake_minimum_required (VERSION 3.10)
project (RawToAces)

set (RawToAces_MAJOR_VERSION 2)
set (RawToAces_MINOR_VERSION 0)
set (RawToAces_PATCH_VERSION 1)
set (RawToAces_VERSION ${RawToAces_MAJOR_VERSION}.${RawToAces_MINOR_VERSION}.${RawToAces_PATCH_VERSION})

set (RawToAcesIdt_LIB "RawToAces_idt")
set (RawToAcesUtil_LIB "RawToAces_util")
set (RawToAcesUtil2_LIB "RawToAces_util2")

set (RawToAces_BIN "rawtoaces")
set (RawToAces2_BIN "rawtoaces2")

set (RawToAces_NAMESPACE "rta")
set (RawToAces_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/src/include")

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    add_compile_options (/W0)
    add_compile_definitions (NOMINMAX)
    add_compile_options ($<$<COMPILE_LANGUAGE:CXX>:/utf-8>)
endif ()

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${warnings}")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${warnings}")
set (CMAKE_DEBUG_POSTFIX "_d" )

# use, i.e. don't skip the full RPATH for the build tree
set (CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set (CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
list (FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if ("${isSystemDir}" STREQUAL "-1")
    set (CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif ("${isSystemDir}" STREQUAL "-1")

# convert install paths to absolute
foreach (p LIB BIN INCLUDE CMAKE)
    set (var INSTALL_${p}_DIR)
    if (NOT IS_ABSOLUTE "${${var}}")
        set (${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
     endif ()
endforeach ()

option (ENABLE_SHARED "Enable Shared Libraries" ON)

if (ENABLE_SHARED)
    set (DO_SHARED SHARED)
else ()
    set (DO_SHARED STATIC)
endif ()

include (configure.cmake)

configure_file(
    src/include/rawtoaces/version.h.in 
    "${PROJECT_BINARY_DIR}/include/rawtoaces/version.h" 
    @ONLY
)
include_directories("${PROJECT_BINARY_DIR}/include")

# generated config files end up in binary dir so to find them, need to add to include path
include_directories (${RawToAces_INCLUDE_DIR})

include (GNUInstallDirs)

add_definitions (
    -DPACKAGE="${PROJECT_NAME}" 
    -DVERSION="${RawToAces_VERSION}"
)
add_subdirectory ("src/rawtoaces_idt")
add_subdirectory ("src/rawtoaces_util")
add_subdirectory ("src/rawtoaces_util2")
add_subdirectory ("src/rawtoaces")
add_subdirectory ("src/rawtoaces2")
 
include(CMakePackageConfigHelpers)

if (PKG_CONFIG_FOUND)
   configure_file (config/RawToAces.pc.in "${PROJECT_BINARY_DIR}/RawToAces.pc" @ONLY)
   install (FILES "${PROJECT_BINARY_DIR}/RawToAces.pc" DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig COMPONENT dev)
endif()

file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake" "# Dummy file")

write_basic_package_version_file(
    "${PROJECT_NAME}ConfigVersion.cmake"
    VERSION ${RawToAces_VERSION}
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    config/${PROJECT_NAME}Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)


install(
    TARGETS ${RawToAcesIdt_LIB} ${RawToAcesUtil_LIB} ${RawToAcesUtil2_LIB} ${RawToAces_BIN} ${RawToAces2_BIN}
    EXPORT ${PROJECT_NAME}Targets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    PUBLIC_HEADER
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
)

install(
    FILES "${PROJECT_BINARY_DIR}/include/rawtoaces/version.h"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/rawtoaces
)

install (
    FILES
        "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake//${PROJECT_NAME}
)

install(
    EXPORT ${PROJECT_NAME}Targets
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
    FILE ${PROJECT_NAME}Targets.cmake
    NAMESPACE ${PROJECT_NAME}::
    EXPORT_LINK_INTERFACE_LIBRARIES
)

install (DIRECTORY data DESTINATION ${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME})

add_subdirectory (unittest)
