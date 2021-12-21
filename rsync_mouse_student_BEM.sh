#!/bin/bash

name='student@158.215.135.160:/home/student'
rsync -vr --exclude "CMake*" --exclude "main" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/obj/* ${name}/researchEISPH/cpp/obj/
rsync -vr --exclude "CMake*" --exclude "main" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/build_sph/* ${name}/researchEISPH/cpp/builds/build_sph/
rsync -vr --exclude "CMake*" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/include/* ${name}/researchEISPH/cpp/include/

# rsync --update -vr --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/*.txt ${name}/research/cpp/builds/
rsync --update -vr --exclude "CMake*" --exclude "main" /Users/tomoaki/Dropbox/markdown/cpp/builds/build_sph/*.obj ${name}/researchEISPH/cpp/builds/build_sph/
