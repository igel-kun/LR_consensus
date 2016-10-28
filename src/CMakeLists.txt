cmake_minimum_required(VERSION 2.6)
project (mast)
# The version number.
set (mast_VERSION_MAJOR 1)
set (mast_VERSION_MINOR 0)

option (DEBUG "debug mode (0=none, 9=max)" OFF)

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

set(CMAKE_COMMON_FLAGS "${CMAKE_COMMON_FLAGS} -pthread -Wall -fmax-errors=2 -ftemplate-backtrace-limit=0 -std=c++11 -DIL_STD")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_COMMON_FLAGS} -O3")

# --- the IMPROVED_FIND_LIBRARY macro (taken from the Bio++ CMakeLists.txt) ---
macro(IMPROVED_FIND_LIBRARY OUTPUT_LIBS lib_name include_to_find)
    FIND_PATH(${lib_name}_INCLUDE_DIR  ${include_to_find} PATHS ENV HOME PATH_SUFFIXES local/bpp-phyl-2.2.0/src local/bpp-seq-2.2.0/src local/bpp-core-2.2.0/src bpp-phyl-2.2.0/src bpp-seq-2.2.0/src bpp-core-2.2.0/src )
    message ( "${include_to_find} ${lib_name} found at ${VAR}" )
  SET(${lib_name}_NAMES ${lib_name} ${lib_name}lib ${lib_name}dll)
  FIND_LIBRARY(${lib_name}_LIBRARY NAMES ${${lib_name}_NAMES} PATHS ENV HOME PATH_SUFFIXES local/bpp-phyl-2.2.0/src local/bpp-seq-2.2.0/src local/bpp-core-2.2.0/src  bpp-phyl-2.2.0/src bpp-seq-2.2.0/src bpp-core-2.2.0/src )

  IF(${lib_name}_LIBRARY)
      MESSAGE("-- Library ${lib_name} found here:")
      MESSAGE("    ${${lib_name}_INCLUDE_DIR}/")
      MESSAGE("    ${${lib_name}_LIBRARY}")
  ELSE(${lib_name}_LIBRARY)
      MESSAGE(FATAL_ERROR "${lib_name} required but not found.")
  ENDIF(${lib_name}_LIBRARY)

  #add the dependency:
  INCLUDE_DIRECTORIES(${${lib_name}_INCLUDE_DIR})
  SET(${OUTPUT_LIBS} ${${OUTPUT_LIBS}} ${${lib_name}_LIBRARY})
ENDMACRO(IMPROVED_FIND_LIBRARY)

# uncomment for static compilation
#set(CMAKE_FIND_LIBRARY_SUFFIXES ".a ${CMAKE_FIND_LIBRARY_SUFFIXES}") 

IMPROVED_FIND_LIBRARY(LIBS bpp-phyl Bpp/Phyl/Tree.h)
IMPROVED_FIND_LIBRARY(LIBS bpp-seq Bpp/Seq/Alphabet/Alphabet.h)
IMPROVED_FIND_LIBRARY(LIBS bpp-core Bpp/Clonable.h)

#Boost libraries:
#find_package(Boost COMPONENTS graph REQUIRED)
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})


message("debug mode (DEBUG)             " ${DEBUG} )
if((${DEBUG} EQUAL 0) OR (${DEBUG} STREQUAL OFF))
  set(CMAKE_BUILD_TYPE RELEASE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNDEBUG")
else((${DEBUG} EQUAL 0) OR (${DEBUG} STREQUAL OFF))
  MESSAGE("building DEBUG version...")
  set(CMAKE_BUILD_TYPE DEBUG)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ddebuglevel=${DEBUG}")
endif((${DEBUG} EQUAL 0) OR (${DEBUG} STREQUAL OFF))



if(STATIC STREQUAL ON)
  SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  SET(BUILD_SHARED_LIBRARIES OFF)
  SET(CMAKE_EXE_LINKER_FLAGS "-static")
endif(STATIC STREQUAL ON)

set(CMAKE_COMMON_FLAGS "${CMAKE_COMMON_FLAGS} -pthread -Wall -fmax-errors=2 -ftemplate-backtrace-limit=0 -std=c++11 -DIL_STD")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_COMMON_FLAGS} -O3 -march=native -Wno-cpp -DNDEBUG -pipe -fno-stack-protector")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_COMMON_FLAGS} -O3 -Wno-cpp -DNDEBUG")


# Source Directory
include_directories("hslinkern")
add_subdirectory ( "hslinkern" )

add_executable( "mast" "mast.cpp" )
target_link_libraries ( "mast" "${LIBS}" "hslinkern" )

add_executable( "mast_all" "mast_all.cpp" )
target_link_libraries ( "mast_all" "${LIBS}" "hslinkern" )

add_executable( "dists" "dists.cpp" )
target_link_libraries ( "dists" "${LIBS}" )



