cmake_minimum_required(VERSION 3.11)

if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(ColourBold "${Esc}[1m")
  set(Red "${Esc}[31m")
  set(Green "${Esc}[32m")
  set(Yellow "${Esc}[33m")
  set(Blue "${Esc}[34m")
  set(Magenta "${Esc}[35m")
  set(Cyan "${Esc}[36m")
  set(White "${Esc}[37m")
  set(BoldRed "${Esc}[1;31m")
  set(BoldGreen "${Esc}[1;32m")
  set(BoldYellow "${Esc}[1;33m")
  set(BoldBlue "${Esc}[1;34m")
  set(BoldMagenta "${Esc}[1;35m")
  set(BoldCyan "${Esc}[1;36m")
  set(BoldWhite "${Esc}[1;37m")
endif()

# Set project name
project(ProjectName)

set(CXX_VERSIONS 13)

# Find C++ compiler
foreach(ver IN LISTS CXX_VERSIONS)
  unset(CXX_COMPILER CACHE)
  string(CONCAT CXX "g++-" ${ver})
  find_program(CXX_COMPILER ${CXX})
  if(CXX_COMPILER)
    message(STATUS "${Green}Found ${CXX}: ${Magenta}${CXX_COMPILER}${ColourReset}")
    set(CMAKE_CXX_COMPILER ${CXX_COMPILER})
    break()
  endif()
endforeach()

if(NOT CXX_COMPILER)
  message(FATAL_ERROR "No suitable g++ version found!")
endif()

# Set common compiler flags
set(CMAKE_CXX_FLAGS "-fopenmp -llapack -lblas -fconcepts")
set(CMAKE_CXX_FLAGS_DEBUG "-Og")
set(CMAKE_CXX_FLAGS_RELEASE "-pthread -Ofast -march=native -funroll-loops -H")

# Set C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")

###########################################################################
# time make
# gold:
# make  41.96s user 1.10s system 100% cpu 42.908 total
# bfd:
# make  43.08s user 1.12s system 100% cpu 44.087 total
# lld"
# make  41.95s user 1.02s system 99% cpu 42.995 total

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fuse-ld=gold")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fuse-ld=gold")

###########################################################################

# Enable ccache
option(CCACHE_ENABLE "If the command ccache is available, use it for compile." ON)
find_program(CCACHE_EXE ccache)

if(CCACHE_ENABLE AND CCACHE_EXE)
  set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_EXE}")
  set(CMAKE_C_COMPILER_LAUNCHER "${CCACHE_EXE}")
  message(STATUS "ccache enabled")
else()
  message(STATUS "ccache not found or disabled")
endif()

set(CMAKE_VERBOSE_MAKEFILE ON)

# Enable parallel build
if(NOT DEFINED CMAKE_BUILD_PARALLEL_LEVEL)
  if(DEFINED CMAKE_HOST_SYSTEM_PROCESSORS AND CMAKE_HOST_SYSTEM_PROCESSORS GREATER 0)
    set(CMAKE_BUILD_PARALLEL_LEVEL ${CMAKE_HOST_SYSTEM_PROCESSORS})
  else()
    # Fallback to platform-specific methods to get processor count
    if(UNIX AND NOT APPLE)
      execute_process(COMMAND nproc OUTPUT_VARIABLE NUM_PROCESSORS OUTPUT_STRIP_TRAILING_WHITESPACE)
    elseif(APPLE)
      execute_process(COMMAND sysctl -n hw.logicalcpu OUTPUT_VARIABLE NUM_PROCESSORS OUTPUT_STRIP_TRAILING_WHITESPACE)
    elseif(WIN32)
      execute_process(COMMAND WMIC CPU Get NumberOfLogicalProcessors OUTPUT_VARIABLE NUM_PROCESSORS OUTPUT_STRIP_TRAILING_WHITESPACE)
    endif()
    set(CMAKE_BUILD_PARALLEL_LEVEL ${NUM_PROCESSORS})
  endif()
  message(STATUS "Parallel build level set to ${CMAKE_BUILD_PARALLEL_LEVEL}")
endif()

# Set default source file
set(DEFAULT_SOURCE_FILE "main.cpp")

if(CMAKE_ARGS)
  list(FIND CMAKE_ARGS "-DSOURCE_FILE" SOURCE_FILE_IDX)

  if(NOT SOURCE_FILE_IDX EQUAL -1)
    list(GET CMAKE_ARGS ${SOURCE_FILE_IDX} SOURCE_FILE_ARG)
    string(REPLACE "-DSOURCE_FILE=" "" SOURCE_FILE ${SOURCE_FILE_ARG})
  endif()
endif()

if(NOT DEFINED SOURCE_FILE)
  set(SOURCE_FILE ${DEFAULT_SOURCE_FILE})
endif()

if(DEFINED OUTPUT_NAME)
  set(BASE_NAME ${OUTPUT_NAME})
else()
  get_filename_component(BASE_NAME ${SOURCE_FILE} NAME_WE) # Extract filename without extension
endif()

message(STATUS "${Green}SOURCE_FILE is ${Magenta}" ${SOURCE_FILE} "${ColourReset}")
message(STATUS "${Green}BASE_NAME is ${Magenta}" ${BASE_NAME} "${ColourReset}")

add_executable(${BASE_NAME} ${SOURCE_FILE})


###########################################################################
# libtet.a のパスを変数に設定
# set(LIBTET_PATH "${CMAKE_CURRENT_SOURCE_DIR}/build_Network/libtet.a")
set(LIBTET_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../include/tetgen1.6.0/libtet.a")
# ファイルの存在を確認
if(EXISTS "${LIBTET_PATH}")
  message(STATUS "Found libtet.a: ${LIBTET_PATH}")
  target_link_libraries(${BASE_NAME} "${LIBTET_PATH}")
else()
  message(WARNING "libtet.a not found at ${LIBTET_PATH}. Skipping linkage.")
endif()
###########################################################################

# Include directories if tetgen.h is not in the standard include path
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Include directories
target_include_directories(${BASE_NAME} PUBLIC
  ${PROJECT_SOURCE_DIR}/../include
  ${PROJECT_SOURCE_DIR}/../../include
)

# Enablce verbose makefile and export compile commands for use with various tools
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


###########################################################################
set(PCH_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/../include/pch.hpp")
set(PCH_OUTPUT "${CMAKE_CURRENT_SOURCE_DIR}/../include/pch.hpp.gch")
add_custom_command(
  OUTPUT ${PCH_OUTPUT}
  COMMAND ${CMAKE_CXX_COMPILER} -std=c++${CMAKE_CXX_STANDARD} -x c++-header -o ${PCH_OUTPUT} ${PCH_HEADER} -Ofast
  DEPENDS ${PCH_HEADER}
  COMMENT "Compiling Precompiled Header: ${PCH_HEADER}"
)
add_custom_target(pch DEPENDS ${PCH_OUTPUT})
add_dependencies(${BASE_NAME} pch)
###########################################################################


message(STATUS "====================================================")

message(STATUS "${Green}CMAKE_CXX_COMPILER is ${Magenta}" ${CMAKE_CXX_COMPILER} "${ColourReset}")
message(STATUS "${Green}CMAKE_CXX_FLAGS is ${Magenta}" ${CMAKE_CXX_FLAGS} "${ColourReset}")
message(STATUS "${Green}CMAKE_CXX_FLAGS_DEBUG is ${Magenta}" ${CMAKE_CXX_FLAGS_DEBUG} "${ColourReset}")
message(STATUS "${Green}CMAKE_CXX_FLAGS_RELEASE is ${Magenta}" ${CMAKE_CXX_FLAGS_RELEASE} "${ColourReset}")
message(STATUS "${Green}CMAKE_CXX_STANDARD is ${Magenta}" ${CMAKE_CXX_STANDARD} "${ColourReset}")
message(STATUS "${Green}CMAKE_CURRENT_SOURCE_DIR is ${Magenta}" ${CMAKE_CURRENT_SOURCE_DIR} "${ColourReset}")
message(STATUS "${Green}CMAKE_CURRENT_BINARY_DIR is ${Magenta}" ${CMAKE_CURRENT_BINARY_DIR} "${ColourReset}")
message(STATUS "${Green}CMAKE_BUILD_PARALLEL_LEVEL is ${Magenta}" ${CMAKE_BUILD_PARALLEL_LEVEL} "${ColourReset}")
message(STATUS "${Green}CMAKE_EXE_LINKER_FLAGS is ${Magenta}" ${CMAKE_EXE_LINKER_FLAGS} "${ColourReset}")
message(STATUS "${Green}CMAKE_SHARED_LINKER_FLAGS is ${Magenta}" ${CMAKE_SHARED_LINKER_FLAGS} "${ColourReset}")
message(STATUS "${Green}CMAKE_CXX_COMPILER_LAUNCHER is ${Magenta}" ${CMAKE_CXX_COMPILER_LAUNCHER} "${ColourReset}")
message(STATUS "${Green}CMAKE_VERBOSE_MAKEFILE is ${Magenta}" ${CMAKE_VERBOSE_MAKEFILE} "${ColourReset}")
message(STATUS "${Green}CMAKE_EXPORT_COMPILE_COMMANDS is ${Magenta}" ${CMAKE_EXPORT_COMPILE_COMMANDS} "${ColourReset}")
message(STATUS "${Green}CMAKE_ARGS is ${Magenta}" ${CMAKE_ARGS} "${ColourReset}")
message(STATUS "${Green}SOURCE_FILE_IDX is ${Magenta}" ${SOURCE_FILE_IDX} "${ColourReset}")
message(STATUS "${Green}SOURCE_FILE_ARG is ${Magenta}" ${SOURCE_FILE_ARG} "${ColourReset}")
message(STATUS "${Green}SOURCE_FILE is ${Magenta}" ${SOURCE_FILE} "${ColourReset}")
message(STATUS "${Green}BASE_NAME is ${Magenta}" ${BASE_NAME} "${ColourReset}")

message(STATUS "${Yellow} ============================================================== ${ColourReset}")
message(STATUS "${Red} MAKE SURE YOU HAVE PRECOMPILED HEADER FILE ``sh precompile.sh`` ${ColourReset}")
message(STATUS "${Yellow} ============================================================== ${ColourReset}")

# /opt/homebrew/bin/g++-13 -std=gnu++2b -x c++-header -fopenmp -fconcepts -pthread -Ofast -march=native -I/Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds -I/Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/../include -I/Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/../../include -o ./pch.hpp.gch ./pch.hpp