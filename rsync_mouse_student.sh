#!/bin/bash

name='student@158.215.135.160:/home/student'
# rsync -vr --exclude "CMake*" --exclude "main" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/obj/* ${name}/research_test/cpp/obj/
# rsync -vr --exclude "CMake*" --exclude "main" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/build_sph/* ${name}/research_test/cpp/builds/build_sph/
rsync -vr --exclude "CMake*" --exclude "main" --exclude "*.vtu" --exclude "*.mov" /Users/tomoaki/Dropbox/markdown/cpp/builds/build_octree/* ${name}/research/cpp/builds/build_octree/
rsync -vr --exclude "CMake*" --exclude "main" --exclude "*.vtu" --exclude "*.mov"  /Users/tomoaki/Dropbox/markdown/cpp/builds/build_lapack/* ${name}/research/cpp/builds/build_lapack/
rsync -vr --exclude "CMake*" --exclude "main" --exclude "*.vtu" --exclude "*.mov"  /Users/tomoaki/Dropbox/markdown/cpp/builds/build_divide_merge/* ${name}/research/cpp/builds/build_divide_merge/
rsync -vr --exclude "CMake*" --exclude "main" --exclude "*.vtu" --exclude "*.mov"  /Users/tomoaki/Dropbox/markdown/cpp/builds/build_bem/* ${name}/research/cpp/builds/build_bem/

rsync -vr --exclude "CMake*" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/include/* ${name}/research/cpp/include/
rsync -vr --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/CMakeLists.txt ${name}/research/cpp/builds/

# rsync --update -vr --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/*.txt ${name}/research/cpp/builds/

# rsync --update -vr --exclude "CMake*" --exclude "main" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/build_dem20210410/* student@158.215.135.160:/home/student/research/cpp/builds/build_dem20210410/
# rsync --update -vr --exclude "CMake*" --exclude "main" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/build_sph/* student@158.215.135.160:/home/student/research/cpp/builds/build_sph/
# rsync --update -vr --exclude "CMake*" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/include/* student@158.215.135.160:/home/student/research/cpp/include/

# rsync --update -vr --exclude "CMake*" --exclude "main" /Users/tomoaki/Dropbox/markdown/cpp/builds/build_sph/*.obj student@158.215.135.160:/home/student/research/cpp/builds/build_sph/
