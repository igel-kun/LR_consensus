cmake_minimum_required(VERSION 2.6)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

set (hslinkern_VERSION_MAJOR 0)
set (hslinkern_VERSION_MINOR 1)
set (hslinkern_NAME \"hslinkern\")

set(CMAKE_CXX_FLAGS "-g3 -Wall -O3")

# add_library(hslinkern STATIC hslinkern.cpp db_map.cpp)
set(SOURCES hslinkern.cpp bitops.cpp db_map.cpp Hypergraph.cpp )
add_library(hslinkern ${SOURCES})

