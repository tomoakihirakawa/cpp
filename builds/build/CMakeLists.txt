# CMakeのバージョンを設定
cmake_minimum_required(VERSION 3.1)
# プロジェクト名と使用する言語を設定
project(test_cmake CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE "-std=c++17 -Ofast -fopenmp")

set(name develop20200527)
set(main ../../source/${name})

# a.outという実行ファイルをmain.cppとhello.cppから作成
add_executable(
  ${name}
  ${main}
  )

target_include_directories(
  ${name}
  PUBLIC
  ${PROJECT_SOURCE_DIR}/../../include
  )

#message(${PROJECT_SOURCE_DIR})
