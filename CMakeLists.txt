cmake_minimum_required(VERSION 3.14)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


project(simple_cellsplit VERSION 1.0 DESCRIPTION "Birth Death with a field."
        LANGUAGES CXX)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")


#include(CheckIPOSupported)
#check_ipo_supported(RESULT supported OUTPUT error)

#if( supported )
#    message(STATUS "IPO / LTO enabled")
#    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
#else()
#	message(STATUS "IPO / LTO not supported: <${error}>")
#endif()


find_package(MPI REQUIRED)
add_definitions(-DOMPI_SKIP_MPICXX)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

find_package(PkgConfig REQUIRED)

set(PKG_CHECK_MODULES_STATIC)

set(CMAKE_FIND_LIBRARY_SUFFIXES .a)


add_subdirectory(src)

MESSAGE(STATUS "Compiler flags:" ${CMAKE_CXX_COMPILE_FLAGS})
MESSAGE(STATUS "Compiler cxx debug flags:" ${CMAKE_CXX_FLAGS_DEBUG})
MESSAGE(STATUS "Compiler cxx release flags:" ${CMAKE_CXX_FLAGS_RELEASE})
MESSAGE(STATUS "Compiler cxx min size flags:" ${CMAKE_CXX_FLAGS_MINSIZEREL})
MESSAGE(STATUS "Compiler cxx flags:" ${CMAKE_CXX_FLAGS})
