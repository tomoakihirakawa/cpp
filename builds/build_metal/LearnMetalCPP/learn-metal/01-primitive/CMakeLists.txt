cmake_minimum_required(VERSION 3.15)
project(MetalApp)

set(CMAKE_CXX_STANDARD 17)

find_library(APPKIT AppKit)
find_library(METAL Metal)
find_library(METALKIT MetalKit)

add_executable(MetalApp 01-primitive.cpp)

target_link_libraries(MetalApp ${APPKIT} ${METAL} ${METALKIT})