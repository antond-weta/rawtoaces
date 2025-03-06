cmake_minimum_required(VERSION 3.10)
include_directories( "${CMAKE_CURRENT_SOURCE_DIR}" )

add_executable (
    ${RawToAces_BIN}
    main.cpp
)

set_property (TARGET ${RawToAces_BIN} PROPERTY CXX_STANDARD 17)

target_link_libraries (
    ${RawToAces_BIN}
    PUBLIC
        ${RawToAcesUtil_LIB}
)

if (LIBRAW_CONFIG_FOUND)
    target_link_libraries (${RawToAces_BIN} PUBLIC libraw::raw)
else ()
    target_link_directories (${RawToAces_BIN} PUBLIC ${libraw_LIBRARY_DIRS})
    target_link_libraries (${RawToAces_BIN} PUBLIC ${libraw_LIBRARIES} ${libraw_LDFLAGS_OTHER})
endif ()

install (TARGETS ${RawToAces_BIN} DESTINATION ${CMAKE_INSTALL_BINDIR} )
