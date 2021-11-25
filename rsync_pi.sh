#!/bin/bash

name='pi@192.168.0.113:/home/pi/'
rsync -vr --exclude "CMakeFiles" --exclude "CMakeCache*" --exclude "main" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/build_pybind11/* ${name}/research/cpp/builds/build_pybind11/
rsync -vr --exclude "CMakeFiles" --exclude "CMakeCache*" --exclude "main" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/include/* ${name}/research/cpp/include/
rsync -vr --exclude "fundamental*" --exclude "CMakeFiles" --exclude "CMakeCache*" --exclude "main" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/python_shared/* ${name}/research/cpp/pi/
