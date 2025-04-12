cmake_minimum_required(VERSION 3.10)
project( RAWTOACES )

set( RTA_MAJOR_VERSION 2 )
set( RTA_MINOR_VERSION 0 )
set( RTA_PATCH_VERSION 0 )
set( RTA_VERSION ${RTA_MAJOR_VERSION}.${RTA_MINOR_VERSION}.${RTA_PATCH_VERSION} )

option( RTA_ENABLE_EIGEN "Enable Eigen dependency" TRUE )
option( RTA_ENABLE_CERES "Enable Ceres-Solver dependency" TRUE )

#set( CMAKE_MACOSX_RPATH 1 )

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
#    set(warnings "/W4 /WX /EHsc")
    add_compile_options ( /W0 )
    add_compile_definitions( NOMINMAX )
endif()

if (NOT CONFIGURED_ONCE)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${warnings}")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${warnings}")
endif()

## Make install directories overridable
set( INSTALL_LIB_DIR lib CACHE PATH "Install directory for libraries" )
set( INSTALL_BIN_DIR bin CACHE PATH "Install directory for executable binaries" )
set( INSTALL_INCLUDE_DIR include CACHE PATH "Install directory for public header files" )
if( WIN32 AND NOT CYGWIN )
  set(DEF_INSTALL_CMAKE_DIR CMake)
else()
  set(DEF_INSTALL_CMAKE_DIR lib/CMake/RAWTOACES)
endif()
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH "Install directory for project CMake files" )

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH FALSE)
# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
# SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
  SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

## convert install paths to absolute
foreach( p LIB BIN INCLUDE CMAKE )
  set( var INSTALL_${p}_DIR )
  if( NOT IS_ABSOLUTE "${${var}}" )
    set( ${var} "${CMAKE_INSTALL_PREFIX}/${${var}}" )
  endif()
endforeach()

option( ENABLE_SHARED "Enable Shared Libraries" ON )
option( RTA_CENTOS7_CERES_HACK "Work around broken config in ceres-solver 1.12" OFF )

if ( ENABLE_SHARED )
  set ( DO_SHARED SHARED )
else ()
  set ( DO_SHARED STATIC )
endif ()

include ( configure.cmake )

add_subdirectory("v1")
add_subdirectory("src")
add_subdirectory("tests")

if ( APPLE OR UNIX )
	install (DIRECTORY data DESTINATION share/rawtoaces)
endif()



