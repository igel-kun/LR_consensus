cmake_minimum_required(VERSION 2.6)
project (mast)
# The version number.
set (mast_VERSION_MAJOR 1)
set (mast_VERSION_MINOR 0)

if(POLICY CMP0012)
  CMAKE_POLICY(SET CMP0012 NEW)
endif(POLICY CMP0012)

set (source_files ${PROJECT_NAME}.cpp AUTHORS BUGS COPYING INSTALL README THANKS)

option (DEBUG "debug mode (0=none, 9=max)" OFF)
option (STATIC "enable static build" OFF)

message("Options: (enable with -D<option>=ON & disable with -D<option>=OFF)")
message("debug mode (DEBUG)             " ${DEBUG} )
message("static build (STATIC)          " ${STATIC} )

add_custom_target ( archive tar -cjf ${PROJECT_NAME}.tar.bz2 ${source_files} )
add_custom_target ( doc doxygen Doxyfile )

add_subdirectory ( "hslinkern" )
add_executable( "${PROJECT_NAME}" "${PROJECT_NAME}.cpp" )

include_directories ("${PROJECT_SOURCE_DIR}" )

# static linking
if(${STATIC} STREQUAL ON)
  set(Boost_USE_STATIC_LIBS   ON)
  SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  SET(BUILD_SHARED_LIBRARIES OFF)
  SET(CMAKE_EXE_LINKER_FLAGS "-static")
endif(${STATIC} STREQUAL ON)

FIND_PACKAGE( Boost COMPONENTS program_options graph REQUIRED )
INCLUDE_DIRECTORIES( "${Boost_INCLUDE_DIR}" )

if(NOT Boost_FOUND)
  MESSAGE(FATAL_ERROR "Could not find Boost.")
endif(NOT Boost_FOUND)

if((${DEBUG} EQUAL 0) OR (${DEBUG} STREQUAL OFF))
  set(CMAKE_BUILD_TYPE RELEASE)
else((${DEBUG} EQUAL 0) OR (${DEBUG} STREQUAL OFF))
  MESSAGE("building DEBUG version...")
  set(CMAKE_BUILD_TYPE DEBUG)
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_COMMON_FLAGS} -Ddebuglevel=${DEBUG}")
endif((${DEBUG} EQUAL 0) OR (${DEBUG} STREQUAL OFF))

set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_COMMON_FLAGS "${CMAKE_COMMON_FLAGS} -pthread -Wall -fmax-errors=2 -ftemplate-backtrace-limit=0 -std=c++11 -DIL_STD")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_COMMON_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_COMMON_FLAGS} -O3 -march=native -DNDEBUG -pipe -m64 -msahf -fno-stack-protector")


