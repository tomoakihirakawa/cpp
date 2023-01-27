#!/bin/bash

# name='student@158.215.135.162:/home/student'
# name='tomoaki@10.0.1.7:~/'
name='student@10.0.1.7:~/'

rsync ~/Dropbox/markdown/cpp/builds/CMakeLists.txt  ${name}/research/cpp/builds/
rsync --update -vr --exclude "*.vtu"  --exclude "*.vtp" --exclude "CMake*" --exclude "main" --exclude "*.mov" --exclude "*.py" ~/Dropbox/markdown/cpp/obj/2022*  ${name}/research/cpp/obj

rsync -vr --exclude "CMake*" --exclude "*.vtu" ~/Dropbox/markdown/cpp/include/* ${name}/research/cpp/include/

rsync --update -vr --exclude "*.vtu" --exclude "*.vtp" --exclude "CMake*" --exclude "main" --exclude "*.mov" --exclude "*.json" ~/Dropbox/markdown/cpp/builds/build_bem/* ${name}/research/cpp/builds/build_bem/

# rsync --update -vr \
# --exclude "*.vtu" --exclude "*.vtp" --exclude "CMake*" --exclude "main" --exclude "*.mov" --exclude "*.json" \
# /Users/tomoaki/Dropbox/markdown/cpp/builds/build_sph/* ${name}/research/cpp/builds/build_sph/
