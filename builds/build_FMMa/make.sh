#! /bin/bash
cp ../../include/fundamental.hpp ./fundamental.cpp
emcc ./fundamental.cpp --bind -s ALLOW_MEMORY_GROWTH=1 -Oz -o embinded.js -v -fopenmp
# emcc ./fundamental.cpp \
#      --bind \
#      --pre-js pre.js\
#      -s TOTAL_MEMORY=512MB\
#      -s WASM=1\
#      -s ALLOW_MEMORY_GROWTH=1 \
#      -Oz \
#      -o embinded.js -v 
