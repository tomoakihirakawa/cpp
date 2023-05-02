# Set minimum CMake version
cmake_minimum_required(VERSION 3.11)

# Set project name
project(ProjectName)

# Check platform
if(APPLE)
  message(STATUS "platform is Apple\n")

  set(CMAKE_CXX_COMPILER "/opt/homebrew/bin/g++-13")

# set(CMAKE_C_COMPILER "/opt/homebrew/bin/gcc-12")
# set(CMAKE_CXX_COMPILER /usr/bin/clang++)
# set(CMAKE_C_COMPILER /usr/bin/clang)
elseif(LINUX)
  message(STATUS "platform is Linux\n")
  set(CMAKE_CXX_COMPILER "/usr/bin/g++-12")

# set(CMAKE_C_COMPILER "/usr/bin/gcc-12")
else()
  message(STATUS "platform is Not Linux\n")
  set(CMAKE_CXX_COMPILER "/usr/bin/g++-12")

  # set(CMAKE_C_COMPILER "/usr/bin/gcc-12")
endif()

# Set common compiler flags
set(CMAKE_CXX_FLAGS "-fopenmp -llapack -lblas -fconcepts -fuse-ld=gold")
set(CMAKE_CXX_FLAGS_DEBUG "-Og")
set(CMAKE_CXX_FLAGS_RELEASE "-pthread -ffast-math -Ofast -Wunused-label -Wunused-value -fno-math-errno -funsafe-math-optimizations -fno-signed-zeros -fno-trapping-math")

# Set C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")

# Check for GNU gold linker
execute_process(
  COMMAND ${CMAKE_C_COMPILER} -fuse-ld=gold -Wl,--version
  ERROR_QUIET
  OUTPUT_VARIABLE ld_version)

if("${ld_version}" MATCHES "GNU gold")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fuse-ld=gold -Wl,--disable-new-dtags")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fuse-ld=gold -Wl,--disable-new-dtags")
endif()

# Set compiler launcher
set(CMAKE_C_COMPILER_LAUNCHER ccache)
set(CMAKE_CXX_COMPILER_LAUNCHER ccache)

# Enable ccache
option(CCACHE_ENABLE "If the command ccache is available, use it for compile." ON)
find_program(CCACHE_EXE ccache)

if(CCACHE_EXE AND CCACHE_ENABLE)
  message(STATUS "ccache is enabled")
  set(CMAKE_C_COMPILER_LAUNCHER "${CMAKE_C_COMPILER_LAUNCHER}" "${CCACHE_EXE}")
  set(CMAKE_CXX_COMPILER_LAUNCHER "${CMAKE_CXX_COMPILER_LAUNCHER}" "${CCACHE_EXE}")
endif()

set(CMAKE_VERBOSE_MAKEFILE ON)

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

# Create executable
message("SOURCE_FILE is " ${SOURCE_FILE})
get_filename_component(BASE_NAME ${SOURCE_FILE} NAME_WE)
add_executable(${BASE_NAME} ${SOURCE_FILE})

# Find LAPACK
find_package(LAPACK)

if(LAPACK_FOUND)
  message(STATUS "LAPACK found!")
  include_directories(${LAPACK_INCLUDE_DIRS})
  target_link_libraries(${BASE_NAME} ${LAPACK_LIBRARIES})
else()
  message(STATUS "LAPACK not found")
endif()

# Set include directories
target_include_directories(${BASE_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/../include)
target_include_directories(${BASE_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/../../include)

message(${PROJECT_SOURCE_DIR})
