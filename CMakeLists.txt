cmake_minimum_required(VERSION 3.1.0)

set (CMAKE_CXX_STANDARD 14)

include(CheckCCompilerFlag)
include(CheckCSourceCompiles)

# Set a default build type if none was specified
set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
        STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
#else()
#    message(STATUS "Build type set to '${CMAKE_BUILD_TYPE}'")
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

option(USE_OPENCL "Use OpenCL support." OFF)
option(BUILD_SSE2_CODEPATHS "(EXPERIMENTAL OPTION, DO NOT DISABLE) Building SSE2-optimized codepaths" ON)

if(BUILD_SSE2_CODEPATHS)
    CHECK_C_COMPILER_FLAG("-msse2" _MSSE2)
    if(NOT _MSSE2)
        MESSAGE(WARNING "Building of SSE2-optimized codepaths is enabled, but the compiler does not understand -msse2.")
        set(BUILD_SSE2_CODEPATHS OFF)
    endif()
endif()

MESSAGE(STATUS "Building SSE2-optimized codepaths: ${BUILD_SSE2_CODEPATHS}")

if (CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
MESSAGE(STATUS "is gnu ${CMAKE_COMPILER_IS_GNUCC} ${CMAKE_COMPILER_IS_GNUCXX}")
MESSAGE(STATUS "cmake c compiler version ${CMAKE_C_COMPILER_VERSION} c++ compiler version ${CMAKE_CXX_COMPILER_VERSION}")
if (${CMAKE_C_COMPILER_VERSION} VERSION_GREATER "6.0.0" OR ${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER "6.0.0")
# Check for __builtin_cpu_supports here due to a bug in clang/llvm. LLVM uses
# libgcc to provide the symbols, so if clang wasn't compiled with libgcc, the
# compiler will say it's available; however the linker cannot find the symbols.
# So do a quick test compile to see if it's actually available.
# See https://llvm.org/bugs/show_bug.cgi?id=25510
check_c_source_compiles("#include <stdio.h>
int main() {
  __builtin_cpu_init();
  printf(\"%d\", __builtin_cpu_supports(\"sse3\"));
}" HAVE_BUILTIN_CPU_SUPPORTS)
if(HAVE_BUILTIN_CPU_SUPPORTS)
    add_definitions("-DHAVE_BUILTIN_CPU_SUPPORTS")
endif(HAVE_BUILTIN_CPU_SUPPORTS)
endif()
MESSAGE(STATUS "Does the compiler support __builtin_cpu_supports(): ${HAVE_BUILTIN_CPU_SUPPORTS}")
endif()

find_package (Threads REQUIRED)
list(APPEND LIBS ${CMAKE_THREAD_LIBS_INIT})

if (NOT OpenCV_DIR)
find_package(OpenCV 3 REQUIRED)
else()
# find package will damage OpenCV_DIR if not find
set(OPENCV_PATH ${OpenCV_DIR})
find_package(OpenCV 4 QUIET PATHS ${OPENCV_PATH} NO_DEFAULT_PATH)
if (NOT OPENCV_FOUND)
find_package(OpenCV 3 REQUIRED PATHS ${OPENCV_PATH} NO_DEFAULT_PATH)
endif(NOT OPENCV_FOUND)
endif(NOT OpenCV_DIR)
if (OPENCV_FOUND)
message(STATUS "Found OpenCV ${OpenCV_VERSION_MAJOR}: ${OpenCV_DIR}")
endif()
include_directories(SYSTEM ${OpenCV_INCLUDE_DIRS})
list(APPEND LIBS ${OpenCV_LIBRARIES})

FIND_PACKAGE(OpenMP)
if (OPENMP_FOUND)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif (OPENMP_FOUND)

find_package(Eigen REQUIRED)
include_directories(SYSTEM ${EIGEN_INCLUDE_DIRS})

if(CMAKE_BUILD_TYPE MATCHES Release)
#IF(MSVC)
#   ADD_DEFINITIONS(/arch:SSE2)
#ENDIF(MSVC)
set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ftree-vectorize -funroll-loops")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ftree-vectorize -funroll-loops")
if (NOT ${CMAKE_CROSSCOMPILING})
set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()
elseif(CMAKE_BUILD_TYPE MATCHES Debug)
message("Debug ${CMAKE_BUILD_TYPE}")
set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g")
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}
                    ${CMAKE_CURRENT_SOURCE_DIR}/util
                    ${CMAKE_CURRENT_SOURCE_DIR}/util/lsd
                    )


add_library(lsd SHARED ${CMAKE_CURRENT_SOURCE_DIR}/util/lsd/lsd.c)

add_library(V SHARED ${CMAKE_CURRENT_SOURCE_DIR}/V.cpp)
target_link_libraries(V lsd)
#target_compile_options(image_reader PRIVATE -Wno-error=deprecated-declarations PRIVATE -Wno-deprecated-declarations)

add_executable(test_v ${CMAKE_CURRENT_SOURCE_DIR}/run.cpp)
target_link_libraries(test_v ${LIBS} V lsd)
