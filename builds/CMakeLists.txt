# ==========================================================================
# CMake Project: BEM_Project for GNU C++ (g++)
#
# 意図: このビルドスクリプトは、Apple Silicon上でHomebrewのg++を使用し、
#       四倍精度演算(__float128)とOpenMPを最大限に活用することを目指します。
#       現代的なCMakeの作法に則り、保守性と移植性を高めています。
# ==========================================================================

cmake_minimum_required(VERSION 3.16)
project(BEM_Project CXX)

# --- ターミナル出力用のカラーコードを設定 ---
if(NOT WIN32 AND CMAKE_COLOR_MAKEFILE)
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(Green "${Esc}[32m")
  set(Magenta "${Esc}[35m")
  set(Yellow "${Esc}[33m")
  set(Red "${Esc}[31m")
  set(BoldCyan "${Esc}[1;36m")
endif()

# ==========================================================================
# コンパイラの指定 (g++-14 を優先)
# ==========================================================================
set(CXX_VERSIONS 15 14 13 12)
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

# ==========================================================================
# 必要なパッケージの探索
# ==========================================================================
# 意図: 必要なライブラリ(OpenMP)をCMakeの公式な方法で探します。
#       これにより、環境ごとのパスの違いなどをCMakeが自動で吸収してくれます。

# find_package(OpenMP REQUIRED)
find_package(OpenCL REQUIRED)

# ==========================================================================
# プロジェクト全体の設定
# ==========================================================================
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_UNITY_BUILD ON) # ユニティビルドを有効化（ビルド高速化のため）
# g++の拡張機能を有効にする (gnu++23相当)
set(CMAKE_CXX_EXTENSIONS ON) 
# コンパイルコマンドのデータベースを生成 (clangdなどのツール連携用)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ccacheを有効化（ビルド高速化のため）
option(CCACHE_ENABLE "Use ccache if available" ON)
if(CCACHE_ENABLE)
    find_program(CCACHE_EXE ccache)
    if(CCACHE_EXE)
        set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_EXE}")
        message(STATUS "ccache enabled")
    endif()
endif()

# ==========================================================================
# 実行ファイルの定義
# ==========================================================================
# 意図: あなたのカスタムワークフローを維持します。
# set(SOURCE_FILE "main.cpp")
# set(BASE_NAME "main")
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

if(DEFINED OUTPUT_NAME)
  set(BASE_NAME ${OUTPUT_NAME})
else()
  get_filename_component(BASE_NAME ${SOURCE_FILE} NAME_WE) # Extract filename without extension
endif()
add_executable(${BASE_NAME} ${SOURCE_FILE})

# ==========================================================================
# ターゲットごとのプロパティ設定
# 意図: CMAKE_CXX_FLAGSのようなグローバル変数ではなく、ターゲットごとに設定を
#       適用するのが現代的なCMakeの作法です。これにより、プロジェクトが
#       複雑になっても設定の管理が容易になります。
# ==========================================================================

# --- インクルードディレクトリ ---
target_include_directories(${BASE_NAME} PUBLIC
  "${PROJECT_SOURCE_DIR}/../include"
  "${PROJECT_SOURCE_DIR}/../../include"
)

# --- コンパイルオプション ---
target_compile_options(${BASE_NAME} PRIVATE
    -fopenmp       # OpenMPを有効化
    -march=native  # Apple Siliconに最適化
    -funroll-loops # ループ展開
    -fconcepts     # C++20 Conceptsを有効化
    -pthread       # マルチスレッド
    # ビルドタイプに応じた最適化レベルをジェネレータ式で指定
    $<$<CONFIG:DEBUG>:-Og>
    $<$<CONFIG:RELEASE>:-O3>
)

# --- ライブラリのリンク ---
# 意図: 必要なライブラリをリンカに指示します。-l<name>はCMakeでは
#       シンプルに名前だけ書くのが作法です。
set(LIBTET_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../include/tetgen1.6.0/libtet.a")
target_link_libraries(${BASE_NAME} PRIVATE
    "${LIBTET_PATH}"
    quadmath     # g++で__float128を使うには必須
    lapack
    blas
    gomp
    OpenCL::OpenCL
    # OpenMP::OpenMP_CXX # find_package(OpenMP)で見つけたライブラリ
)

# ==========================================================================
# 最終メッセージ出力
# ==========================================================================
# 意図: ビルド設定の要約を分かりやすく表示し、次のアクションを促します。
string(REPEAT "=" 70 DIVIDER)
message(STATUS "")
message(STATUS "${Yellow}${DIVIDER}${ColourReset}")
message(STATUS "${BoldCyan}  Project Configuration Summary for g++${ColourReset}")
message(STATUS "${Yellow}${DIVIDER}${ColourReset}")
message(STATUS "${Green}Target Name       : ${Magenta}${BASE_NAME}${ColourReset}")
message(STATUS "${Green}C++ Compiler      : ${Magenta}${CMAKE_CXX_COMPILER}${ColourReset}")
message(STATUS "${Green}C++ Standard      : ${Magenta}gnu++${CMAKE_CXX_STANDARD}${ColourReset}")
message(STATUS "${Green}Build Type        : ${Magenta}${CMAKE_BUILD_TYPE}${ColourReset}")
message(STATUS "${Red}  MAKE SURE YOU HAVE COMPILED tetgen1.6.0 BEFORE BUILDING!${ColourReset}")
message(STATUS "${Yellow}To build, run: ${BoldCyan}cmake --build . -j$(sysctl -n hw.logicalcpu)${ColourReset}")
message(STATUS "${Yellow}${DIVIDER}${ColourReset}")
# cmake_minimum_required(VERSION 3.16)

# if(NOT WIN32)
#   string(ASCII 27 Esc)
#   set(ColourReset "${Esc}[m")
#   set(ColourBold "${Esc}[1m")
#   set(Red "${Esc}[31m")
#   set(Green "${Esc}[32m")
#   set(Yellow "${Esc}[33m")
#   set(Blue "${Esc}[34m")
#   set(Magenta "${Esc}[35m")
#   set(Cyan "${Esc}[36m")
#   set(White "${Esc}[37m")
#   set(BoldRed "${Esc}[1;31m")
#   set(BoldGreen "${Esc}[1;32m")
#   set(BoldYellow "${Esc}[1;33m")
#   set(BoldBlue "${Esc}[1;34m")
#   set(BoldMagenta "${Esc}[1;35m")
#   set(BoldCyan "${Esc}[1;36m")
#   set(BoldWhite "${Esc}[1;37m")
# endif()

# # Set project name
# project(ProjectName)

# set(CXX_VERSIONS 14)

# ###########################################################################

# # Find C++ compiler
# # foreach(ver IN LISTS CXX_VERSIONS)
# #   unset(CXX_COMPILER CACHE)
# #   string(CONCAT CXX "g++-" ${ver})
# #   find_program(CXX_COMPILER ${CXX})
# #   if(CXX_COMPILER)
# #     message(STATUS "${Green}Found ${CXX}: ${Magenta}${CXX_COMPILER}${ColourReset}")
# #     set(CMAKE_CXX_COMPILER ${CXX_COMPILER})
# #     break()
# #   endif()
# # endforeach()
# # if(NOT CXX_COMPILER)
# #   message(FATAL_ERROR "No suitable g++ version found!")
# # endif()

# # HomebrewでインストールしたLLVMのclang++を探す
# find_program(CXX_COMPILER clang++ PATHS "/opt/homebrew/opt/llvm/bin" NO_DEFAULT_PATH)
# if(CXX_COMPILER)
#     message(STATUS "${Green}Found Homebrew clang++: ${Magenta}${CXX_COMPILER}${ColourReset}")
#     set(CMAKE_CXX_COMPILER ${CXX_COMPILER})
#     # OpenMPライブラリへのパスを追加（リンカエラーが出る場合）
#     find_library(OMP_LIB omp PATHS "/opt/homebrew/opt/libomp/lib")
#     if(OMP_LIB)
#         message(STATUS "${Green}Found OpenMP library: ${Magenta}${OMP_LIB}${ColourReset}")
#         # これは直接リンクするより、-Lフラグでパスを示す方が良い場合が多い
#         # target_link_libraries(${BASE_NAME} PRIVATE ${OMP_LIB})
#         link_directories("/opt/homebrew/opt/libomp/lib")
#         include_directories("/opt/homebrew/opt/libomp/include")
#     endif()
# else()
#     # 見つからない場合はシステムのデフォルトclang++を使用（OpenMPは使えない可能性）
#     find_program(CXX_COMPILER clang++)
#     message(FATAL_ERROR "${Yellow}Homebrew clang++ not found. Please run 'brew install llvm'. Using system default is not recommended.${ColourReset}")
#     set(CMAKE_CXX_COMPILER ${CXX_COMPILER})
# endif()

# ##########################################################################

# # Set common compiler flags
# set(CMAKE_CXX_FLAGS "-fopenmp -llapack -lblas -fconcepts")
# set(CMAKE_CXX_FLAGS_DEBUG "-pthread -Og -march=native -funroll-loops -lquadmath")
# set(CMAKE_CXX_FLAGS_RELEASE "-pthread -O3 -march=native -funroll-loops -lquadmath")

# # Set C++ standard
# set(CMAKE_CXX_STANDARD 23)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")

# ###########################################################################
# # time make
# # gold:
# # make  41.96s user 1.10s system 100% cpu 42.908 total
# # bfd:
# # make  43.08s user 1.12s system 100% cpu 44.087 total
# # lld"
# # make  41.95s user 1.02s system 99% cpu 42.995 total
# # set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
# # set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}")
# ###########################################################################

# # Enable ccache
# option(CCACHE_ENABLE "If the command ccache is available, use it for compile." ON)
# find_program(CCACHE_EXE ccache)

# if(CCACHE_ENABLE AND CCACHE_EXE)
#   set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_EXE}")
#   set(CMAKE_C_COMPILER_LAUNCHER "${CCACHE_EXE}")
#   message(STATUS "ccache enabled")
# else()
#   message(STATUS "ccache not found or disabled")
# endif()

# set(CMAKE_VERBOSE_MAKEFILE ON)

# # Enable parallel build
# if(NOT DEFINED CMAKE_BUILD_PARALLEL_LEVEL)
#   if(DEFINED CMAKE_HOST_SYSTEM_PROCESSORS AND CMAKE_HOST_SYSTEM_PROCESSORS GREATER 0)
#     set(CMAKE_BUILD_PARALLEL_LEVEL ${CMAKE_HOST_SYSTEM_PROCESSORS})
#   else()
#     # Fallback to platform-specific methods to get processor count
#     if(UNIX AND NOT APPLE)
#       execute_process(COMMAND nproc OUTPUT_VARIABLE NUM_PROCESSORS OUTPUT_STRIP_TRAILING_WHITESPACE)
#     elseif(APPLE)
#       execute_process(COMMAND sysctl -n hw.logicalcpu OUTPUT_VARIABLE NUM_PROCESSORS OUTPUT_STRIP_TRAILING_WHITESPACE)
#     elseif(WIN32)
#       execute_process(COMMAND WMIC CPU Get NumberOfLogicalProcessors OUTPUT_VARIABLE NUM_PROCESSORS OUTPUT_STRIP_TRAILING_WHITESPACE)
#     endif()
#     set(CMAKE_BUILD_PARALLEL_LEVEL ${NUM_PROCESSORS})
#   endif()
#   message(STATUS "Parallel build level set to ${CMAKE_BUILD_PARALLEL_LEVEL}")
# endif()

# # Set default source file
# set(DEFAULT_SOURCE_FILE "main.cpp")

# if(CMAKE_ARGS)
#   list(FIND CMAKE_ARGS "-DSOURCE_FILE" SOURCE_FILE_IDX)

#   if(NOT SOURCE_FILE_IDX EQUAL -1)
#     list(GET CMAKE_ARGS ${SOURCE_FILE_IDX} SOURCE_FILE_ARG)
#     string(REPLACE "-DSOURCE_FILE=" "" SOURCE_FILE ${SOURCE_FILE_ARG})
#   endif()
# endif()

# if(NOT DEFINED SOURCE_FILE)
#   set(SOURCE_FILE ${DEFAULT_SOURCE_FILE})
# endif()

# if(DEFINED OUTPUT_NAME)
#   set(BASE_NAME ${OUTPUT_NAME})
# else()
#   get_filename_component(BASE_NAME ${SOURCE_FILE} NAME_WE) # Extract filename without extension
# endif()

# message(STATUS "${Green}SOURCE_FILE is ${Magenta}" ${SOURCE_FILE} "${ColourReset}")
# message(STATUS "${Green}BASE_NAME is ${Magenta}" ${BASE_NAME} "${ColourReset}")

# add_executable(${BASE_NAME} ${SOURCE_FILE})


# ###########################################################################
# # libtet.a のパスを変数に設定
# # set(LIBTET_PATH "${CMAKE_CURRENT_SOURCE_DIR}/build_Network/libtet.a")
# set(LIBTET_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../include/tetgen1.6.0/libtet.a")
# # ファイルの存在を確認
# if(EXISTS "${LIBTET_PATH}")
#   message(STATUS "Found libtet.a: ${LIBTET_PATH}")
#   target_link_libraries(${BASE_NAME} "${LIBTET_PATH}")
# else()
#   message(WARNING "libtet.a not found at ${LIBTET_PATH}. Skipping linkage.")
# endif()
# ###########################################################################

# # Include directories if tetgen.h is not in the standard include path
# include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# # Include directories
# target_include_directories(${BASE_NAME} PUBLIC
#   ${PROJECT_SOURCE_DIR}/../include
#   ${PROJECT_SOURCE_DIR}/../../include
# )

# # Enablce verbose makefile and export compile commands for use with various tools
# set(CMAKE_VERBOSE_MAKEFILE ON)
# set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


# ###########################################################################
# # set(PCH_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/../include/pch.hpp")
# # set(PCH_OUTPUT "${CMAKE_CURRENT_SOURCE_DIR}/../include/pch.hpp.gch")
# # add_custom_command(
# #   OUTPUT ${PCH_OUTPUT}
# #   COMMAND ${CMAKE_CXX_COMPILER} -std=c++${CMAKE_CXX_STANDARD} -x c++-header -o ${PCH_OUTPUT} ${PCH_HEADER} -Ofast
# #   DEPENDS ${PCH_HEADER}
# #   COMMENT "Compiling Precompiled Header: ${PCH_HEADER}"
# # )
# # add_custom_target(pch DEPENDS ${PCH_OUTPUT})
# # add_dependencies(${BASE_NAME} pch)
# ###########################################################################


# message(STATUS "====================================================")

# message(STATUS "${Green}CMAKE_CXX_COMPILER is ${Magenta}" ${CMAKE_CXX_COMPILER} "${ColourReset}")
# message(STATUS "${Green}CMAKE_CXX_FLAGS is ${Magenta}" ${CMAKE_CXX_FLAGS} "${ColourReset}")
# message(STATUS "${Green}CMAKE_CXX_FLAGS_DEBUG is ${Magenta}" ${CMAKE_CXX_FLAGS_DEBUG} "${ColourReset}")
# message(STATUS "${Green}CMAKE_CXX_FLAGS_RELEASE is ${Magenta}" ${CMAKE_CXX_FLAGS_RELEASE} "${ColourReset}")
# message(STATUS "${Green}CMAKE_CXX_STANDARD is ${Magenta}" ${CMAKE_CXX_STANDARD} "${ColourReset}")
# message(STATUS "${Green}CMAKE_CURRENT_SOURCE_DIR is ${Magenta}" ${CMAKE_CURRENT_SOURCE_DIR} "${ColourReset}")
# message(STATUS "${Green}CMAKE_CURRENT_BINARY_DIR is ${Magenta}" ${CMAKE_CURRENT_BINARY_DIR} "${ColourReset}")
# message(STATUS "${Green}CMAKE_BUILD_PARALLEL_LEVEL is ${Magenta}" ${CMAKE_BUILD_PARALLEL_LEVEL} "${ColourReset}")
# message(STATUS "${Green}CMAKE_EXE_LINKER_FLAGS is ${Magenta}" ${CMAKE_EXE_LINKER_FLAGS} "${ColourReset}")
# message(STATUS "${Green}CMAKE_SHARED_LINKER_FLAGS is ${Magenta}" ${CMAKE_SHARED_LINKER_FLAGS} "${ColourReset}")
# message(STATUS "${Green}CMAKE_CXX_COMPILER_LAUNCHER is ${Magenta}" ${CMAKE_CXX_COMPILER_LAUNCHER} "${ColourReset}")
# message(STATUS "${Green}CMAKE_VERBOSE_MAKEFILE is ${Magenta}" ${CMAKE_VERBOSE_MAKEFILE} "${ColourReset}")
# message(STATUS "${Green}CMAKE_EXPORT_COMPILE_COMMANDS is ${Magenta}" ${CMAKE_EXPORT_COMPILE_COMMANDS} "${ColourReset}")
# message(STATUS "${Green}CMAKE_ARGS is ${Magenta}" ${CMAKE_ARGS} "${ColourReset}")
# message(STATUS "${Green}SOURCE_FILE_IDX is ${Magenta}" ${SOURCE_FILE_IDX} "${ColourReset}")
# message(STATUS "${Green}SOURCE_FILE_ARG is ${Magenta}" ${SOURCE_FILE_ARG} "${ColourReset}")
# message(STATUS "${Green}SOURCE_FILE is ${Magenta}" ${SOURCE_FILE} "${ColourReset}")
# message(STATUS "${Green}BASE_NAME is ${Magenta}" ${BASE_NAME} "${ColourReset}")

# message(STATUS "${Yellow} ============================================================== ${ColourReset}")
# message(STATUS "${Red} MAKE SURE YOU HAVE PRECOMPILED HEADER FILE ``sh precompile.sh`` ${ColourReset}")
# message(STATUS "${Yellow} ============================================================== ${ColourReset}")

# /opt/homebrew/bin/g++-13 -std=gnu++2b -x c++-header -fopenmp -fconcepts -pthread -Ofast -march=native -I/Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds -I/Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/../include -I/Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/../../include -o ./pch.hpp.gch ./pch.hpp