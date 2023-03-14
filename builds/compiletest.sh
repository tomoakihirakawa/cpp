#!/bin/bash
RED="\033[1;31m"
NC='\033[0m'
DIR='/Users/tomoaki/Dropbox/code/cpp/builds'

echo "============================================"
tmp=${DIR}"/build_computaional_geometry"
(
    echo ${RED}"check compile"${NC}
    echo ${RED}${tmp}${NC}
    cd ${tmp}
    make
)
echo "============================================"
tmp=${DIR}"/build_sph"
(
    echo ${RED}"check compile"${NC}
    echo ${RED}${tmp}${NC}
    cd ${tmp}
    make
)
echo "============================================"
tmp=${DIR}"/build_octree"
(
    echo ${RED}"check compile"${NC}
    echo ${RED}${tmp}${NC}
    cd ${tmp}
    make
)
echo "============================================"
tmp=${DIR}"/build_octree"
(
    echo ${RED}"check compile"${NC}
    echo ${RED}${tmp}${NC}
    cd ${tmp}
    make
)
echo "============================================"
tmp=${DIR}"/build_divide_merge"
(
    echo ${RED}"check compile"${NC}
    echo ${RED}${tmp}${NC}
    cd ${tmp}
    make
)
echo "============================================"
echo "done"