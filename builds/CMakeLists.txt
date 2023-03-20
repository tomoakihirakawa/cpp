# CMakeのバージョンを設定
cmake_minimum_required(VERSION 3.11)

# プロジェクト名と使用する言語を設定
# project(network_bem CXX)
project(ProjectName)

# message(${CMAKE_SOURCE_DIR}) set (CMAKE_MODULE_PATH
# "${CMAKE_SOURCE_DIR}/CMake") include(cotire)

# ##############################################################################
if(APPLE)
  message(STATUS ">> Apple")

  set(CMAKE_CXX_COMPILER "/opt/homebrew/bin/g++-12")
  set(CMAKE_CXX_FLAGS "-fopenmp -llapack -lblas -fconcepts -fuse-ld=gold")
  set(CMAKE_CXX_FLAGS_DEBUG "-Og")

  # set(CMAKE_CXX_FLAGS "-Ofast -fopenmp") set(CMAKE_CXX_FLAGS " -llapack -lblas
  # -lgfortran -include ${PROJECT_SOURCE_DIR}/../include/fundamental.hpp
  # -std=c++17 -fuse-ld=gold" )

  # set(CMAKE_CXX_COMPILER "/usr/bin/clang++")
  # set(CMAKE_CXX_FLAGS " -fopenmp=libomp -llapack -lblas -fuse-ld=gold")
  # set(CMAKE_CXX_FLAGS_DEBUG "-Og")

  # set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -fopenmp -Wnonnull -Wnonnull-compare -Wopenmp-simd -Wparentheses -Wpointer-sign    -Wreorder     -Wreturn-type    -Wsequence-point   -Wsign-compare  -Wsizeof-pointer-memaccess   -Wstrict-aliasing -Wstrict-overflow=1 -Wswitch -Wtautological-compare -Wtrigraphs  -Wuninitialized    -Wunknown-pragmas  -Wunused-function  -Wunused-label      -Wunused-value -Wunused-variable -Wvolatile-registercm-var")
  # set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -fopenacc -fopenmp -Wunused-function -Wunused-label -Wunused-value -Wunused-variable")
  set(CMAKE_CXX_FLAGS_RELEASE "-pthread -ffast-math -Ofast -Wunused-label -Wunused-value -fno-math-errno -funsafe-math-optimizations -fno-signed-zeros -fno-trapping-math")

# ############################################################################
elseif(LINUX)
  set(CMAKE_CXX_COMPILER "/usr/bin/g++-12")
  message(STATUS ">> Linux")
  set(CMAKE_CXX_FLAGS " -fopenmp -llapack -lblas  -fconcepts -fuse-ld=gold")
  set(CMAKE_CXX_FLAGS_DEBUG "-Og")
  set(CMAKE_CXX_FLAGS_RELEASE "-pthread -ffast-math -Ofast -Wunused-label -Wunused-value -fno-math-errno -funsafe-math-optimizations -fno-signed-zeros -fno-trapping-math")

# ############################################################################
else()
  set(CMAKE_CXX_COMPILER "/usr/bin/g++-12")
  message(STATUS ">> Not Linux")
  set(CMAKE_CXX_FLAGS "-fopenmp  -llapack -lblas -fconcepts -fuse-ld=gold")
  set(CMAKE_CXX_FLAGS_DEBUG "-Og")
  set(CMAKE_CXX_FLAGS_RELEASE "-pthread -ffast-math -Ofast -Wunused-label -Wunused-value -fno-math-errno -funsafe-math-optimizations -fno-signed-zeros -fno-trapping-math")
endif()

# ##############################################################################
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")

#
execute_process(
  COMMAND ${CMAKE_C_COMPILER} -fuse-ld=gold -Wl,--version
  ERROR_QUIET
  OUTPUT_VARIABLE ld_version)

if("${ld_version}" MATCHES "GNU gold")
  set(CMAKE_EXE_LINKER_FLAGS
    "${CMAKE_EXE_LINKER_FLAGS} -fuse-ld=gold -Wl,--disable-new-dtags")
  set(CMAKE_SHARED_LINKER_FLAGS
    "${CMAKE_SHARED_LINKER_FLAGS} -fuse-ld=gold -Wl,--disable-new-dtags")
endif()

set(CMAKE_C_COMPILER_LAUNCHER ccache)
set(CMAKE_CXX_COMPILER_LAUNCHER ccache)

option(CCACHE_ENABLE "If the command ccache is avilable, use it for compile." ON)

find_program(CCACHE_EXE ccache)

if(CCACHE_EXE)
  if(CCACHE_ENABLE)
    message(STATUS ">> Enable ccache")

    if(CMAKE_C_COMPILER_LAUNCHER)
      set(CMAKE_C_COMPILER_LAUNCHER "${CMAKE_C_COMPILER_LAUNCHER}"
        "${CCACHE_EXE}")
    else()
      set(CMAKE_C_COMPILER_LAUNCHER "${CCACHE_EXE}")
    endif()

    if(CMAKE_CXX_COMPILER_LAUNCHER)
      set(CMAKE_CXX_COMPILER_LAUNCHER "${CMAKE_CXX_COMPILER_LAUNCHER}"
        "${CCACHE_EXE}")
    else()
      set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_EXE}")
    endif()
  endif()
endif()

set(CMAKE_VERBOSE_MAKEFILE ON)

# set(main main.cpp) set(out main)
# mainという実行ファイルをmain.cppから作成
# set(DEFAULT_MAIN_FILE "main.cpp")
# if(NOT DEFINED MAIN_FILE)
# set(MAIN_FILE ${DEFAULT_MAIN_FILE})
# endif()
# add_executable(main ${MAIN_FILE})
#
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

message(${SOURCE_FILE})
get_filename_component(BASE_NAME ${SOURCE_FILE} NAME_WE)
add_executable(${BASE_NAME} ${SOURCE_FILE})

# cotire(main)
# if(APPLE)

# else() # * ------------------------------------------------------ */
# find_package(BLAS) find_package(LAPACK) if(LAPACK_FOUND AND BLAS_FOUND)
# set(lapackblas_libraries "${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}") endif()
# target_link_libraries(main PUBLIC ${lapackblas_libraries}) # *
# ------------------------------------------------------ */ endif()
target_include_directories(${BASE_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/../include)
target_include_directories(${BASE_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/../../include)

message(${PROJECT_SOURCE_DIR})

# memo
# https://edgarluque.com/blog/cmake-precompiled-headersを見る限りMyProjectがはじめにくる
