#!/bin/zsh

#@ Network.hppに関するヘッダーファイルは，このプリコンパイル`pch.hpp`に含めていない
#@ 変更があまりされないものだけを，`pch.hpp`に含めている

SCRIPT_DIR=$(cd $(dirname $0) && pwd)

# (
#    cd "${SCRIPT_DIR}/tetgen1.6.0"
#    sh clean
#    cmake -DCMAKE_BUILD_TYPE=Release ./
#    make
# )

(
   cd "${SCRIPT_DIR}"
   /opt/homebrew/bin/g++-13 -std=gnu++2b -x c++-header -fopenmp -fconcepts -pthread -Ofast -march=native -I/Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds -I/Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/../include -I/Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/../../include -o ./pch.hpp.gch ./pch.hpp
   # /opt/homebrew/bin/g++-13 -v -Wall -std=gnu++2b -x c++-header -fopenmp -fconcepts -pthread -Ofast -march=native -I"${SCRIPT_DIR}" -I"${SCRIPT_DIR}/tetgen" -o ./pch.hpp.gch ./pch.hpp
)
