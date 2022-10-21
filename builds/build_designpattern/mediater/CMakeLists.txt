# CMakeのバージョンを設定
cmake_minimum_required(VERSION 3.1)
# プロジェクト名と使用する言語を設定
project(test_cmake CXX)
set(CMAKE_CXX_COMPILER /usr/local/bin/g++-9)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE "-std=c++17 -Wall -Og")

set(main pattern.cpp)
set(out pattern)

# a.outという実行ファイルをmain.cppとhello.cppから作成
add_executable(
  ${out}
  ${main}
  )

target_include_directories(
  ${out}
  PUBLIC
  ${PROJECT_SOURCE_DIR}/../../../include
  )

message(${PROJECT_SOURCE_DIR})
