#!/bin/bash

name='student@158.215.135.161:/home/student'
# rsync -vr --exclude "CMake*" --exclude "main" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/obj/* ${name}/researchEISPH/cpp/obj/
rsync -vr --exclude "CMake*" --exclude "main" --exclude "*.vtu" --exclude "*.mov" --exclude "*.json" /Users/tomoaki/Dropbox/markdown/cpp/builds/build_bem/* ${name}/research/cpp/builds/build_bem/

rsync -vr --exclude "CMake*" --exclude "main" --exclude "*.vtu" --exclude "*.mov" /Users/tomoaki/Dropbox/markdown/cpp/include/* ${name}/research/cpp/include/
# rsync -vr --exclude "CMake*" --exclude "main" --exclude "*.vtu" --exclude "*.mov" /Users/tomoaki/Dropbox/markdown/cpp/obj/2022* ${name}/research/cpp/obj/

# rsync --update -vr --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/*.txt ${name}/research/cpp/builds/
# rsync --update -vr --exclude "CMake*" --exclude "main" /Users/tomoaki/Dropbox/markdown/cpp/builds/build_sph/*.obj ${name}/researchEISPH/cpp/builds/build_sph/
