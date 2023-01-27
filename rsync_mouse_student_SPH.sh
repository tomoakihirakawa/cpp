#!/bin/bash

name='student@158.215.135.160:/home/student'
# rsync -vr --exclude "CMake*" --exclude "main" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/obj/* ${name}/research/cpp/obj/
# rsync -vr --exclude "CMake*" --exclude "main" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/build_bem/* ${name}/research/cpp/builds/build_bem/
# rsync -vr --exclude "CMake*" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/include/* ${name}/research/cpp/include/

# rsync --update -vr --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/*.txt ${name}/research/cpp/builds/
# rsync --update -vr --exclude "*.vtu" --exclude "main" --exclude "*.mov"/Users/tomoaki/Dropbox/markdown/cpp/CMake* ${name}/research/cpp/
rsync -vr --exclude --exclude "*.vtu" --exclude "*.vtp" --exclude "*.vtu" "CMake*" --exclude "main" --exclude "*.mov" /Users/tomoaki/Dropbox/markdown/cpp/include/* ${name}/research/cpp/include/
rsync --update -vr --exclude "*.vtu" --exclude "*.vtp" --exclude "CMake*" --exclude "main" --exclude "*.mov" /Users/tomoaki/Dropbox/markdown/cpp/builds/build_sph/* ${name}/research/cpp/builds/build_sph/
# rsync --update -vr --exclude "CMake*" --exclude "main" /Users/tomoaki/Dropbox/markdown/cpp/builds/build_sph/*.obj ${name}/research/cpp/builds/build_sph/
rsync --update -vr --exclude "CMake*" --exclude "main" /Users/tomoaki/Dropbox/markdown/cpp/obj/2023* ${name}/research/cpp/obj

# scp student@158.215.135.160:/home/student/research/cpp/builds/build_octree/OC5* ./
