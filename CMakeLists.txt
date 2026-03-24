#SPDX-License-Identifier: MIT
#Copyright (c) Gavin W. Wilson

cmake_minimum_required (VERSION 3.30)
project (scdepth C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(PYBIND11_FINDPYTHON ON)
set(CMAKE_CXX_SCAN_FOR_MODULES OFF)


SET( EXTERNAL_DIR "${PROJECT_SOURCE_DIR}/extern" )
include(ExternalProject)

add_subdirectory(${EXTERNAL_DIR}/pybind11)

MESSAGE("${EXTERNAL_DIR}" )
ExternalProject_Add(htslib
  SOURCE_DIR        ${EXTERNAL_DIR}/htslib
  PREFIX            ${EXTERNAL_DIR}/htslib-prefix
  URL               https://github.com/samtools/htslib/releases/download/1.22.1/htslib-1.22.1.tar.bz2
  URL_HASH          SHA256=3dfa6eeb71db719907fe3ef7c72cb2ec9965b20b58036547c858c89b58c342f7
  CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env
                    CC=${CMAKE_C_COMPILER}
                    CFLAGS=-g\ -Wall\ -O2\ -fvisibility=hidden\ -fPIC
                    ./configure --disable-libcurl
  BUILD_COMMAND     ${MAKE_EXE}
  BUILD_IN_SOURCE   1
  INSTALL_COMMAND   ${CMAKE_COMMAND} -E echo "Skipping install step."
  BUILD_BYPRODUCTS  ${EXTERNAL_DIR}/htslib/libhts.a
)

ExternalProject_Get_Property(htslib source_dir)
add_library(htslib_imp STATIC IMPORTED GLOBAL)
set_target_properties(htslib_imp PROPERTIES
  IMPORTED_LOCATION             "${source_dir}/libhts.a"
  INTERFACE_INCLUDE_DIRECTORIES "${source_dir}"
)
add_dependencies(htslib_imp htslib)

MESSAGE( "BUILD TYPE: ${CMAKE_BUILD_TYPE}" )
MESSAGE( "INSTALL PREFIX: ${CMAKE_INSTALL_PREFIX}" )

# Build the link list
find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)
find_package(BZip2 QUIET)
find_package(LibLZMA QUIET)
find_package(ZSTD QUIET)
find_library(LIBDEFLATE_LIBRARY NAMES deflate libdeflate QUIET)

set(_hts_deps ZLIB::ZLIB)
if(TARGET ZSTD::ZSTD)
    MESSAGE("Found ZSTD Lib")
    list(APPEND _hts_deps ZSTD::ZSTD)
endif()
if(BZip2_FOUND)
    MESSAGE("Found BZip2 Lib")
    list(APPEND _hts_deps BZip2::BZip2)
endif()
if(LibLZMA_FOUND)
    MESSAGE("Found LZMA Lib")
    list(APPEND _hts_deps LibLZMA::LibLZMA)
endif()
if(LIBDEFLATE_LIBRARY)
    MESSAGE("Found deflate Lib")
    list(APPEND _hts_deps ${LIBDEFLATE_LIBRARY})
endif()

set_property(TARGET htslib_imp APPEND PROPERTY
  INTERFACE_LINK_LIBRARIES "${_hts_deps}"
)

#MESSAGE("HTS LIBS = ${_hts_deps}")

add_subdirectory("src")
