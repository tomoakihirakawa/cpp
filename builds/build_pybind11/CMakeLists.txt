cmake_minimum_required(VERSION 3.11)

if(APPLE)
  set(CMAKE_CXX_COMPILER /opt/homebrew/bin/g++-13)
endif()

# プロジェクト名と使用する言語を設定
project(my_pybind11 CXX)

set(CMAKE_CXX_FLAGS "-Ofast -std=c++2b")
set(CMAKE_CXX_STANDARD 23)

if(APPLE)
  message(STATUS ">>> Apple")
  add_subdirectory(pybind11)
  pybind11_add_module(LightHillRobot_pybind SHARED LightHillRobot_pybind.cpp)
elseif(LINUX)
  message(STATUS ">>> Linux")
  add_subdirectory(pybind11)
  pybind11_add_module(LightHillRobot_pybind SHARED LightHillRobot_pybind.cpp)
else()
  message(STATUS ">>> Not Linux")
  add_subdirectory(pybind11)
  pybind11_add_module(LightHillRobot_pybind SHARED LightHillRobot_pybind.cpp)
endif()
