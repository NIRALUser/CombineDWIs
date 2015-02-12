CMAKE_MINIMUM_REQUIRED( VERSION 3.0 )
CMAKE_POLICY(VERSION 3.0)
project(CombineDWIs)

find_package(SlicerExecutionModel REQUIRED)
include( ${SlicerExecutionModel_USE_FILE} )
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})


set(Filter CombineDWIs)

#snippet taken from: http://www.guyrutenberg.com/2014/01/05/enabling-c11-c0x-in-cmake/
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
  set( EXTRA_FLAG "-std=c++11" )
elseif(COMPILER_SUPPORTS_CXX0X)
  set( EXTRA_FLAG "-std=c++0x" )
else()
   message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()
if( EXTRA_FLAG )
  set( CMAKE_CXX_FLAGS_LIST ${CMAKE_CXX_FLAGS} )
  separate_arguments(CMAKE_CXX_FLAGS_LIST)
  get_property( helpstring CACHE CMAKE_CXX_FLAGS PROPERTY HELPSTRING )
  list(FIND CMAKE_CXX_FLAGS_LIST ${EXTRA_FLAG} found_flag )
  if( found_flag LESS 0 )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_FLAG}" CACHE STRING "${helpstring}" FORCE )
  endif()
endif()
SEMMacroBuildCLI(
    NAME ${Filter}
    EXECUTABLE_ONLY
    TARGET_LIBRARIES ${ITK_LIBRARIES}
)
