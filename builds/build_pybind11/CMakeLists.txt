cmake_minimum_required(VERSION 3.11)

if(APPLE)
  set(CMAKE_CXX_COMPILER /opt/homebrew/bin/g++-13)
endif()

# プロジェクト名と使用する言語を設定
project(my_pybind11 CXX)

set(CMAKE_CXX_FLAGS "-Ofast -std=c++2b")
set(CMAKE_CXX_STANDARD 23)

# Set default input and output file
set(INPUT "main.cpp" CACHE STRING "Input file")
set(OUTPUT "main" CACHE STRING "Output file")

if(APPLE)
  message(STATUS ">>> Apple")
  add_subdirectory(pybind11)
  pybind11_add_module(${OUTPUT} SHARED ${INPUT})
elseif(LINUX)
  message(STATUS ">>> Linux")
  add_subdirectory(pybind11)
  pybind11_add_module(${OUTPUT} SHARED ${INPUT})
else()
  message(STATUS ">>> Not Linux")
  add_subdirectory(pybind11)
  pybind11_add_module(${OUTPUT} SHARED ${INPUT})
endif()
