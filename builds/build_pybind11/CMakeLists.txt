cmake_minimum_required(VERSION 3.11)

# cmake -DCMAKE_BUILD_TYPE=Release ./ -DINPUT=LighthillRobot.cpp -DOUTPUT=LighthillRobot -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-13

# プロジェクト名と使用する言語を設定
project(my_pybind11 CXX)

set(CMAKE_CXX_FLAGS "-Ofast -std=c++2b")

# Set default input and output file
set(INPUT "main.cpp" CACHE STRING "Input file")
set(OUTPUT "main" CACHE STRING "Output file")

if(APPLE)
    message(STATUS ">>> Apple")
    set(CMAKE_CXX_COMPILER /opt/homebrew/bin/g++-13)
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
