#!/bin/bash
RED="\033[1;31m"
NC='\033[0m'

tmp="/Users/tomoaki/Dropbox/markdown/cpp/builds/build_geometry"
(
    echo ${RED}check compile${NC}
    echo ${RED}${tmp}${NC}
    cd ${tmp}
    make
)

tmp="/Users/tomoaki/Dropbox/markdown/cpp/builds/build_sph"
(
    echo ${RED}check compile${NC}
    echo ${RED}${tmp}${NC}
    cd ${tmp}
    make
)

tmp="/Users/tomoaki/Dropbox/markdown/cpp/builds/build_octree"
(
    echo ${RED}check compile${NC}
    echo ${RED}${tmp}${NC}
    cd ${tmp}
    make
)

tmp="/Users/tomoaki/Dropbox/markdown/cpp/builds/build_geometry"
(
    echo ${RED}check compile${NC}
    echo ${RED}${tmp}${NC}
    cd ${tmp}
    make
)

tmp="/Users/tomoaki/Dropbox/markdown/cpp/builds/build_divide_merge"
(
    echo ${RED}check compile${NC}
    echo ${RED}${tmp}${NC}
    cd ${tmp}
    make
)

echo "done"
