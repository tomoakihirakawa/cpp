# CMakeのバージョンを設定
cmake_minimum_required(VERSION 3.11)

# プロジェクト名と使用する言語を設定
project(network_bem CXX)

# message(${CMAKE_SOURCE_DIR}) set (CMAKE_MODULE_PATH
# "${CMAKE_SOURCE_DIR}/CMake") include(cotire)

# ##############################################################################
if(APPLE)
  set(CMAKE_CXX_COMPILER g++-12)
  message(STATUS ">> Apple")

  # set(CMAKE_CXX_FLAGS "-Ofast -fopenmp") set(CMAKE_CXX_FLAGS " -llapack -lblas
  # -lgfortran -include ${PROJECT_SOURCE_DIR}/../include/fundamental.hpp
  # -std=c++17 -fuse-ld=gold" )
  set(CMAKE_CXX_FLAGS " -llapack -lblas -fconcepts -fuse-ld=gold")
  set(CMAKE_CXX_FLAGS_DEBUG "-Og")

  # set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -fopenmp -Wnonnull -Wnonnull-compare -Wopenmp-simd -Wparentheses -Wpointer-sign    -Wreorder     -Wreturn-type    -Wsequence-point   -Wsign-compare  -Wsizeof-pointer-memaccess   -Wstrict-aliasing -Wstrict-overflow=1 -Wswitch -Wtautological-compare -Wtrigraphs  -Wuninitialized    -Wunknown-pragmas  -Wunused-function  -Wunused-label      -Wunused-value -Wunused-variable -Wvolatile-register-var")
  # set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -fopenmp -Wunused-function -Wunused-label -Wunused-value -Wunused-variable")
  set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -fopenmp -Wunused-label -Wunused-value")

# ############################################################################
elseif(LINUX)
  set(CMAKE_CXX_COMPILER "/usr/bin/g++-11")
  message(STATUS ">> Linux")
  set(CMAKE_CXX_FLAGS " -llapack -lblas -fuse-ld=gold")
  set(CMAKE_CXX_FLAGS_DEBUG "-Og")
  set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -fopenmp -march=native")

# ############################################################################
else()
  set(CMAKE_CXX_COMPILER "/usr/bin/g++-11")
  message(STATUS ">> Not Linux")
  set(CMAKE_CXX_FLAGS " -llapack -lblas -fuse-ld=gold")
  set(CMAKE_CXX_FLAGS_DEBUG "-Og")
  set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -fopenmp -march=native")
endif()

# ##############################################################################

# if(UNIX OR APPLE)
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

# endif()
set(CMAKE_C_COMPILER_LAUNCHER ccache)
set(CMAKE_CXX_COMPILER_LAUNCHER ccache)
set(CMAKE_CXX_STANDARD 20)

option(CCACHE_ENABLE "If the command ccache is avilable, use it for compile."
  ON)

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
add_executable(main main.cpp)

# cotire(main)

# if(APPLE)

# else() # * ------------------------------------------------------ */
# find_package(BLAS) find_package(LAPACK) if(LAPACK_FOUND AND BLAS_FOUND)
# set(lapackblas_libraries "${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}") endif()
# target_link_libraries(main PUBLIC ${lapackblas_libraries}) # *
# ------------------------------------------------------ */ endif()
target_include_directories(main PUBLIC ${PROJECT_SOURCE_DIR}/../include)
target_include_directories(main PUBLIC ${PROJECT_SOURCE_DIR}/../../include)

message(${PROJECT_SOURCE_DIR})

# memo
# https://edgarluque.com/blog/cmake-precompiled-headersを見る限りMyProjectがはじめにくる
