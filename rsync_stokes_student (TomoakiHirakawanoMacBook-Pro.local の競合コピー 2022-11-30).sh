#!/bin/bash

# name='student@158.215.135.161:/home/student'
name='student@10.0.1.8:/home/student'

rsync --update -vr --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/*.txt ${name}/research/cpp/builds/
rsync -vr --exclude --exclude "*.vtu" --exclude "*.vtp" --exclude "*.vtu" "CMake*" --exclude "main" --exclude "*.mov" /Users/tomoaki/Dropbox/markdown/cpp/include/* ${name}/research/cpp/include/
rsync --update -vr --exclude "*.vtu" --exclude "*.vtp" --exclude "CMake*" --exclude "main" --exclude "*.mov" /Users/tomoaki/Dropbox/markdown/cpp/builds/build_sph/* ${name}/research/cpp/builds/build_sph/
rsync --update -vr --exclude "*.vtu" --exclude "*.vtp" --exclude "CMake*" --exclude "main" --exclude "*.mov" /Users/tomoaki/Dropbox/markdown/cpp/builds/build_bem/* ${name}/research/cpp/builds/build_bem/
rsync --update -vr --exclude "*.vtu" --exclude "*.vtp" --exclude "CMake*" --exclude "main" --exclude "*.mov" /Users/tomoaki/Dropbox/markdown/cpp/builds/build_divide_merge/* ${name}/research/cpp/builds/build_divide_merge/
rsync --update -vr --exclude "CMake*" --exclude "main" /Users/tomoaki/Dropbox/markdown/cpp/obj/2022*/* ${name}/research/cpp/obj/
rsync -vr --exclude "CMake*" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/include/* ${name}/research/cpp/include/
# rsync -vr --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/CMakeLists.txt ${name}/research/cpp/builds/

# rsync --update -vr --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/*.txt ${name}/research/cpp/builds/

# rsync --update -vr --exclude "CMake*" --exclude "main" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/build_dem20210410/* student@158.215.135.160:/home/student/research/cpp/builds/build_dem20210410/
# rsync --update -vr --exclude "CMake*" --exclude "main" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/build_sph/* student@158.215.135.160:/home/student/research/cpp/builds/build_sph/
# rsync --update -vr --exclude "CMake*" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/include/* student@158.215.135.160:/home/student/research/cpp/include/

# rsync --update -vr --exclude "CMake*" --exclude "main" /Users/tomoaki/Dropbox/markdown/cpp/builds/build_sph/*.obj student@158.215.135.160:/home/student/research/cpp/builds/build_sph/
