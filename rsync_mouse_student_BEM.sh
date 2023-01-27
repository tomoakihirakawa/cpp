#!/bin/bash
name='student@158.215.135.160:/home/student'

# rsync -vr --exclude "CMake*" --exclude "main" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/obj/* ${name}/researchEISPH/cpp/obj/
rsync -vr --exclude "CMake*" --exclude "*mov" --exclude "main" --exclude "*.vtp" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/build_bem/* ${name}/research/cpp/builds/build_bem/
rsync -vr --exclude "CMake*" --exclude "*mov" --exclude "main" --exclude "*.vtp" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/build_divide_merge/* ${name}/research/cpp/builds/build_divide_merge/
rsync -vr --exclude "CMake*" --exclude "*mov" --exclude "*.vtp" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/include/* ${name}/research/cpp/include/
# rsync -vr --exclude "CMake*" --exclude "*mov" --exclude "*.vtp" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/obj/2022* ${name}/research/cpp/obj/
# rsync --update -vr --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/*.txt ${name}/research/cpp/builds/
# rsync --update -vr --exclude "CMake*" --exclude "main" /Users/tomoaki/Dropbox/markdown/cpp/builds/build_sph/*.obj ${name}/researchEISPH/cpp/builds/build_sph/
