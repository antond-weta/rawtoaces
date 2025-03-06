# Copyright Contributors to the rawtoaces project.
# SPDX-License-Identifier: Apache-2.0
# https://github.com/AcademySoftwareFoundation/rawtoaces

include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")

add_executable (${RawToAces2_BIN}
    main.cpp
)

set_property (TARGET ${RawToAces2_BIN} PROPERTY CXX_STANDARD 17)

target_include_directories (
    ${RawToAces2_BIN}
    PUBLIC
        ${OpenImageIO_INCLUDE_DIR}
)

target_link_libraries (
    ${RawToAces2_BIN}
    PUBLIC
        ${RawToAcesUtil2_LIB}
)

set_target_properties (${RawToAces2_BIN} PROPERTIES XCODE_SCHEME_ARGUMENTS "${CMAKE_XCODE_SCHEME_ARGUMENTS}")

install (TARGETS ${RawToAces2_BIN} DESTINATION bin)

