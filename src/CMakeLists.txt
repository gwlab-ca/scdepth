#SPDX-License-Identifier: MIT
#Copyright (c) Gavin W. Wilson

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src/)

add_compile_options(-Wall -Wextra -Wpedantic)

set(SOURCES
    tags.cpp
    gzstream.cpp
    shards.cpp
    counter.cpp
    downsampler.cpp
    dsdata.cpp
    dsthread.cpp
    umis.cpp
    goverlaps.cpp
    gtf.cpp
    visbins.cpp
)

set(BINDING_SOURCES
    bindings.cpp
    bind_ds.cpp
    bind_vis.cpp
    bind_cntr.cpp
    bind_govs.cpp
)

#add_executable(scdepth-bin scdepth.cpp ${SOURCES})
#add_dependencies(scdepth-bin htslib)
#target_link_libraries(scdepth-bin PUBLIC htslib_imp  Threads::Threads ZLIB::ZLIB)
#target_include_directories(scdepth-bin SYSTEM PUBLIC ${ZLIB_INCLUDE_DIRS})
#set_target_properties(scdepth-bin PROPERTIES OUTPUT_NAME scdepth)

pybind11_add_module(_bindings ${BINDING_SOURCES} ${SOURCES})
target_link_libraries(_bindings PUBLIC htslib_imp  Threads::Threads ZLIB::ZLIB)
target_include_directories(_bindings SYSTEM PUBLIC ${ZLIB_INCLUDE_DIRS})
add_dependencies(_bindings htslib)
install(TARGETS _bindings LIBRARY DESTINATION scdepth/bindings)
