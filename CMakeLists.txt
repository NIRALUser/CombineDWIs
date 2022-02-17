# ==== Define cmake build policies that affect compilation and linkage default behaviors
#
# Set the NEWEST_VALIDATED_POLICIES_VERSION string to the newest cmake version
# policies that provide successful builds. By setting NEWEST_VALIDATED_POLICIES_VERSION
# to a value greater than the oldest policies, all policies between
# OLDEST_VALIDATED_POLICIES_VERSION and CMAKE_VERSION (used for this build)
# are set to their NEW behaivor, thereby suppressing policy warnings related to policies
# between the OLDEST_VALIDATED_POLICIES_VERSION and CMAKE_VERSION.
#
# CMake versions greater than the NEWEST_VALIDATED_POLICIES_VERSION policies will
# continue to generate policy warnings "CMake Warning (dev)...Policy CMP0XXX is not set:"
#
set(OLDEST_VALIDATED_POLICIES_VERSION "3.16.3")
set(NEWEST_VALIDATED_POLICIES_VERSION "3.19.7")
cmake_minimum_required(VERSION ${OLDEST_VALIDATED_POLICIES_VERSION}...${NEWEST_VALIDATED_POLICIES_VERSION} FATAL_ERROR)

#-----------------------------------------------------------------------------
# Enable C++14
#-----------------------------------------------------------------------------
if(NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 14)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
message(STATUS "${_msg} - C++${CMAKE_CXX_STANDARD}")
if(NOT CMAKE_CXX_STANDARD MATCHES "^(14|17|20)$")
  message(FATAL_ERROR "CMAKE_CXX_STANDARD must be set to 14, 17, or 20")
endif()

project(CombineDWIs)

find_package(SlicerExecutionModel REQUIRED)
include( ${SlicerExecutionModel_USE_FILE} )
find_package(ITK 5.3 REQUIRED)
include(${ITK_USE_FILE})


set(Filter CombineDWIs)

#snippet taken from: http://www.guyrutenberg.com/2014/01/05/enabling-c11-c0x-in-cmake/
include(CheckCXXCompilerFlag)
if( EXTRA_FLAG )
  set( CMAKE_CXX_FLAGS_LIST ${CMAKE_CXX_FLAGS} )
  separate_arguments(CMAKE_CXX_FLAGS_LIST)
  get_property( helpstring CACHE CMAKE_CXX_FLAGS PROPERTY HELPSTRING )
  list(FIND CMAKE_CXX_FLAGS_LIST ${EXTRA_FLAG} found_flag )
  if( found_flag LESS 0 )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_FLAG}" CACHE STRING "${helpstring}" FORCE )
  endif()
endif()

if(NOT COMMAND SETIFEMPTY)
  macro(SETIFEMPTY)
    set(KEY ${ARGV0})
    set(VALUE ${ARGV1})
    if(NOT ${KEY})
      set(${ARGV})
    endif()
  endmacro()
endif()

SETIFEMPTY(INSTALL_RUNTIME_DESTINATION bin)
SETIFEMPTY(INSTALL_LIBRARY_DESTINATION bin)
SETIFEMPTY(INSTALL_ARCHIVE_DESTINATION lib/static)
SETIFEMPTY(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
SETIFEMPTY(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SETIFEMPTY(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)

SEMMacroBuildCLI(
  NAME ${Filter}
  TARGET_LIBRARIES ${ITK_LIBRARIES}
  EXECUTABLE_ONLY
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
  ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}
  INSTALL_RUNTIME_DESTINATION ${INSTALL_RUNTIME_DESTINATION}
  INSTALL_LIBRARY_DESTINATION ${INSTALL_LIBRARY_DESTINATION}
  INSTALL_ARCHIVE_DESTINATION ${INSTALL_ARCHIVE_DESTINATION}
)
