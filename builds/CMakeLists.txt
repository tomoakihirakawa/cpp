# Set minimum CMake version
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

set(CXX_VERSIONS 15 14 13 12)

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
set(CMAKE_CXX_FLAGS "-fopenmp -llapack -lblas -fconcepts -fuse-ld=gold -funroll-loops")
set(CMAKE_CXX_FLAGS_DEBUG "-Og")
set(CMAKE_CXX_FLAGS_RELEASE "-pthread -Ofast -Wunused-label -Wunused-value -fno-math-errno -funsafe-math-optimizations -fno-signed-zeros -fno-trapping-math -march=native -ftree-vectorize -fopt-info-vec-optimized  -funroll-loops")


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

# # Use LLVM's lld linker
# if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
#   set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fuse-ld=lld")
#   set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fuse-ld=lld")
# endif()

# Set compiler launcher
set(CMAKE_C_COMPILER_LAUNCHER ccache)
set(CMAKE_CXX_COMPILER_LAUNCHER ccache)

# Enable ccache
option(CCACHE_ENABLE "If the command ccache is available, use it for compile." ON)
find_program(CCACHE_EXE ccache)

# if(CCACHE_EXE AND CCACHE_ENABLE)
#   message(STATUS "${Green}ccache is enabled${ColourReset}")
#   set(CMAKE_C_COMPILER_LAUNCHER "${CMAKE_C_COMPILER_LAUNCHER}" "${CCACHE_EXE}")
#   set(CMAKE_CXX_COMPILER_LAUNCHER "${CMAKE_CXX_COMPILER_LAUNCHER}" "${CCACHE_EXE}")
# endif()

if(CCACHE_ENABLE AND CCACHE_EXE)
  set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_EXE}")
  set(CMAKE_C_COMPILER_LAUNCHER "${CCACHE_EXE}")
  message(STATUS "ccache enabled")
else()
  message(STATUS "ccache not found or disabled")
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
# message(STATUS "${Green}SOURCE_FILE is ${Magenta}" ${SOURCE_FILE} "${ColourReset}")
# get_filename_component(BASE_NAME ${SOURCE_FILE} NAME_WE)
# add_executable(${BASE_NAME} ${SOURCE_FILE})

# # Set default output file name
# set(DEFAULT_OUTPUT_NAME "MyExecutable")

# Allow the output file name to be passed as a command-line option
if(DEFINED OUTPUT_NAME)
  set(BASE_NAME ${OUTPUT_NAME})
else()
  get_filename_component(BASE_NAME ${SOURCE_FILE} NAME_WE) # Extract filename without extension
endif()

message(STATUS "${Green}SOURCE_FILE is ${Magenta}" ${SOURCE_FILE} "${ColourReset}")
message(STATUS "${Green}BASE_NAME is ${Magenta}" ${BASE_NAME} "${ColourReset}")

add_executable(${BASE_NAME} ${SOURCE_FILE})

# Find and link LAPACK BRARIES})

# Directly specify the path to libtet.a for linking
target_link_libraries(${BASE_NAME} "${CMAKE_CURRENT_SOURCE_DIR}/build_Network/libtet.a")

# Include directories if tetgen.h is not in the standard include path
include_directories(${CMAKE_CURRENT_SOURCE_DIR})


# # Locate and link Eigen3 if find inlcudes and libraries only if found
# find_package(Eigen3 REQUIRED NO_MODULE)
# target_link_libraries(${BASE_NAME} PUBLIC Eigen3::Eigen)

# # Find and link CGAL
# find_package(CGAL REQUIRED)
# target_link_libraries(${BASE_NAME} PUBLIC CGAL::CGAL)

# find_package(CGAL QUIET) # Try to find CGAL quietly

# if (CGAL_FOUND)
#     message(STATUS "CGAL found")
#     include(${CGAL_USE_FILE})
#     target_link_libraries(${BASE_NAME} ${CGAL_LIBRARIES} gmp mpfr)
#     set(CMAKE_CXX_FLAGS "-I/opt/homebrew/Cellar/gmp/6.3.0/include -I/path/to/CGAL/include -L/opt/homebrew/Cellar/gmp/6.3.0/lib -lgmp -lCGAL -lCGAL_Core")
# else()
#     message(STATUS "CGAL not found, proceeding without it")
#     # Link with other necessary libraries (if any), but not with CGAL
#     # You might still need to link with GMP and MPFR if they are used independently of CGAL
#     # target_link_libraries(${BASE_NAME} gmp mpfr)
# endif()

# Include directories
target_include_directories(${BASE_NAME} PUBLIC
  ${PROJECT_SOURCE_DIR}/../include
  ${PROJECT_SOURCE_DIR}/../../include
)


# Enablce verbose makefile and export compile commands for use with various tools
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)