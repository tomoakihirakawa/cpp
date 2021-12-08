# CMakeのバージョンを設定

cmake_minimum_required(VERSION 3.11)
# プロジェクト名と使用する言語を設定
project(my_pybind11 CXX)

# set(CMAKE_CXX_COMPILER /usr/local/bin/g++-10)

# message(${CMAKE_SOURCE_DIR}) set (CMAKE_MODULE_PATH
# "${CMAKE_SOURCE_DIR}/CMake") include(cotire)

if(APPLE)
  # set(CMAKE_CXX_COMPILER g++-10)
  message(STATUS ">>> Apple")
  set(CMAKE_CXX_FLAGS "-Ofast")
  # set(CMAKE_CXX_FLAGS "-g0") set(CMAKE_CXX_FLAGS "-include
  # ${PROJECT_SOURCE_DIR}/../include/fundamental.hpp") set(CMAKE_CXX_FLAGS_DEBUG
  # "-Og -H") set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -fopenmp")

  set(-DCMAKE_CXX_STANDARD=17)
  # set(CMAKE_CXX_FLAGS "-Wall -Wextra -Ofast")
  set(CMAKE_CXX_FLAGS "-Ofast -std=c++17")
  add_subdirectory(pybind11)
  pybind11_add_module(fundamental SHARED fundamental.cpp)
elseif(LINUX)
  # set(CMAKE_CXX_COMPILER /usr/bin/g++)
  message(STATUS ">>> Linux")
  set(PYBIND11_CPP_STANDARD -std=c++17)
  # set(CMAKE_CXX_FLAGS "-Wall -Wextra -Og")
  add_subdirectory(pybind11)
  pybind11_add_module(fundamental SHARED fundamental.cpp)
else()
  # デフォルトのコンパイラをg++にしたい場合 export CXX=/usr/bin/g++-8を実行しておく
  # set(CMAKE_CXX_COMPILER g++-8)
  message(STATUS ">>> Not Linux")
  set(-DCMAKE_CXX_STANDARD=17)
  set(CMAKE_CXX_FLAGS "-Ofast -std=c++17")
  add_subdirectory(pybind11)

  pybind11_add_module(fundamental SHARED fundamental.cpp)

  # メモ
  #
  # sudo apt install python3.9-dev
  #
  # ln -s /usr/lib/python3.9 /usr/include/python3.9
  #
  # git clone https://github.com/pybind/pybind11.git
  #
  # リアルタイムプロットがうまくいかない，使っているpythonバージョンを変えてみる
endif()
