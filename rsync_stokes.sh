#!/bin/bash

name='stokes@10.0.1.2:/home/stokes'
rsync --update -vr --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/*.txt ${name}/research/cpp/builds/
# rsync --update -vr --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/*.txt ${name}/research/cpp/builds/build_sph/
rsync --update -vr --exclude "CMake*" --exclude "main" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/obj/* ${name}/research/cpp/obj/
rsync --update -vr --exclude "CMake*" --exclude "main" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/build_bem/* ${name}/research/cpp/builds/build_bem/
rsync --update -vr --exclude "CMake*" --exclude "main" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/build_lapack/* ${name}/research/cpp/builds/build_lapack/
rsync --update -vr --exclude "CMake*" --exclude "main" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/build_sph/* ${name}/research/cpp/builds/build_sph/
rsync --update -vr --exclude "CMakeC*" --exclude "main" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/build_pybind11/* ${name}/research/cpp/builds/build_pybind11/
rsync --update -vr --exclude "CMakeC*" --exclude "main" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/python_shared/* ${name}/research/python_shared/
rsync --update -vr --exclude "CMake*" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/include/* ${name}/research/cpp/include/


# #!/bin/bash

# # # name='stokes@158.215.135.161:/home/stokes'
# # name='stokes@10.0.1.10:/home/stokes'
# # rsync --update -vr --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/*.txt ${name}/research/cpp/builds/
# # # rsync --update -vr --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/*.txt ${name}/research/cpp/builds/build_sph/
# # rsync --update -vr --exclude "CMake*" --exclude "main" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/obj/* ${name}/research/cpp/obj/
# # rsync --update -vr --exclude "CMake*" --exclude "main" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/build_bem/* ${name}/research/cpp/builds/build_bem/
# # rsync --update -vr --exclude "CMake*" --exclude "main" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/build_lapack/* ${name}/research/cpp/builds/build_lapack/
# # rsync --update -vr --exclude "CMake*" --exclude "main" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/build_sph/* ${name}/research/cpp/builds/build_sph/
# # rsync --update -vr --exclude "CMakeC*" --exclude "main" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/build_pybind11/* ${name}/research/cpp/builds/build_pybind11/
# # rsync --update -vr --exclude "CMakeC*" --exclude "main" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/python_shared/* ${name}/research/python_shared/
# # rsync --update -vr --exclude "CMake*" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/include/* ${name}/research/cpp/include/


# name='stokes@10.0.1.10:/home/stokes'
# rsync --update -vr --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/*.txt ${name}/research/cpp/builds/
# # rsync --update -vr --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/*.txt ${name}/research/cpp/builds/build_sph/
# rsync --update -vr --exclude "CMake*" --exclude "main" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/obj/* ${name}/research/cpp/obj/
# rsync --update -vr --exclude "CMake*" --exclude "main" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/build_bem/* ${name}/research/cpp/builds/build_bem/
# rsync --update -vr --exclude "CMake*" --exclude "main" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/build_lapack/* ${name}/research/cpp/builds/build_lapack/
# # rsync --update -vr --exclude "CMake*" --exclude "main" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/build_sph/* ${name}/research/cpp/builds/build_sph/
# # rsync --update -vr --exclude "CMakeC*" --exclude "main" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/builds/build_pybind11/* ${name}/research/cpp/builds/build_pybind11/
# # rsync --update -vr --exclude "CMakeC*" --exclude "main" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/python_shared/* ${name}/research/python_shared/
# rsync --update -vr --exclude "CMake*" --exclude "*.vtu" /Users/tomoaki/Dropbox/markdown/cpp/include/* ${name}/research/cpp/include/
