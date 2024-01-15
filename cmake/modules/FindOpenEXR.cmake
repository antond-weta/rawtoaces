#
# A simple cmake find module for openexr
#

find_package(PkgConfig QUIET)
if(PKG_CONFIG_FOUND)
  pkg_check_modules(PC_OPENEXR QUIET openexr)
endif()

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  # Under Mac OS, if the user has installed using brew, the package config
  # information is not entirely accurate in that it refers to the
  # true install path but brew maintains links to /usr/local for you
  # which they want you to use
  if(PC_OPENEXR_FOUND AND "${PC_OPENEXR_INCLUDEDIR}" MATCHES ".*/Cellar/.*")
    set(_openexr_HINT_INCLUDE /usr/local/include)
    set(_openexr_HINT_LIB /usr/local/lib)
  endif()
endif()

if(PC_OPENEXR_FOUND)
  set(openexr_CFLAGS ${PC_OPENEXR_CFLAGS_OTHER})
  set(openexr_LIBRARY_DIRS ${PC_OPENEXR_LIBRARY_DIRS})
  set(openexr_LDFLAGS ${PC_OPENEXR_LDFLAGS_OTHER})
  if("${_openexr_HINT_INCLUDE}" STREQUAL "")
    set(_openexr_HINT_INCLUDE ${PC_OPENEXR_INCLUDEDIR} ${PC_OPENEXR_INCLUDE_DIRS})
    set(_openexr_HINT_LIB ${PC_OPENEXR_LIBDIR} ${PC_OPENEXR_LIBRARY_DIRS})
  endif()
else()
  if(UNIX)
    set(openexr_CFLAGS "-pthread")
    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    else()
      set(openexr_LDFLAGS "-pthread")
    endif()
  endif()
endif()

find_path(openexr_INCLUDE_DIR ImfVersion.h
          HINTS ${_openexr_HINT_INCLUDE}
          PATH_SUFFIXES openexr )
if(openexr_INCLUDE_DIR AND EXISTS "${openexr_INCLUDE_DIR}/openexr_version.h")
    message("openexr_INCLUDE_DIR ${openexr_INCLUDE_DIR}")
    
    
    set(openexr_VERSION ${PC_OPENEXR_VERSION})
    
        message("openexr_VERSION ${openexr_VERSION}")

    if("${openexr_VERSION}" STREQUAL "")
      file(STRINGS "${openexr_INCLUDE_DIR}/openexr_version.h" openexr_version_str
           REGEX "^#define[\t ]+OPENEXR_VERSION_STR[\t ]+\".*")

      string(REGEX REPLACE "^#define[\t ]+OPENEXR_VERSION[\t ]+\"([^ \\n]*)\".*"
             "\\1" openexr_VERSION "${openexr_version_str}")
      # message (STATUS "${openexr_version_str}")
      unset(openexr_version_str)
    endif()
endif()

find_library(EXR_LIBRARY
             NAMES OpenEXR
             HINTS ${_openexr_HINT_LIB}
)

if(EXR_LIBRARY)

        message("EXR_LIBRARY ${EXR_LIBRARY}")
        
  set(openexr_LIBRARY ${EXR_LIBRARY})
  mark_as_advanced(${EXR_LIBRARY})
endif()

unset(_openexr_HINT_INCLUDE)
unset(_openexr_HINT_LIB)

set( openexr_LIBRARIES ${openexr_LIBRARY} )
set( openexr_INCLUDE_DIRS ${openexr_INCLUDE_DIR} )

if(NOT PC_OPENEXR_FOUND)
get_filename_component(openexr_LDFLAGS_OTHER ${openexr_LIBRARY} PATH)
set(openexr_LDFLAGS_OTHER -L${openexr_LDFLAGS_OTHER})
endif()

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set OpenEXR_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(OpenEXR
                                  REQUIRED_VARS openexr_LIBRARY openexr_INCLUDE_DIR
                                  VERSION_VAR openexr_VERSION
                                  FAIL_MESSAGE "Unable to find openexr library" )

# older versions of cmake don't support FOUND_VAR to find_package_handle
# so just do it the hard way...
if(OPENEXR_FOUND AND NOT openexr_FOUND)
  set(openexr_FOUND 1)
endif()

mark_as_advanced(openexr_INCLUDE_DIR openexr_LIBRARY )
