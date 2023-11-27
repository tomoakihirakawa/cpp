#! /bin/bash

sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=remesh.cpp
make

# ./remesh ~/Dropbox/code/cpp/obj/Ren2015/wavemaker.obj ~/Dropbox/code/cpp/obj/Ren2015 wavemaker 100
# ./remesh ~/Dropbox/code/cpp/obj/Ren2015/float.obj ~/Dropbox/code/cpp/obj/Ren2015 float 100
# ./remesh ~/Dropbox/code/cpp/obj/Ren2015/water.obj ~/Dropbox/code/cpp/obj/Ren2015 water 400
# ./remesh ~/Dropbox/code/cpp/obj/Ren2015/water.obj ~/Dropbox/code/cpp/obj/Ren2015 water 100
# ./remesh ~/Dropbox/code/cpp/obj/Ren2015/tank.obj ~/Dropbox/code/cpp/obj/Ren2015 tank 100

# define string variable

# ---------------------------------------------------------------------------- #
# path='/Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/obj/Ren2015_multiple/'
# ./remesh ${path}float_a.obj ${path} float_a 100
# ./remesh ${path}float_b.obj ${path} float_b 100
# ./remesh ${path}float_c.obj ${path} float_c 100
# ./remesh ${path}float_d.obj ${path} float_d 100
# ./remesh ${path}float_e.obj ${path} float_e 100
# ./remesh ${path}float_f.obj ${path} float_f 100
# ./remesh ${path}float_g.obj ${path} float_g 100
# ./remesh ${path}wavemaker.obj ${path} wavemaker 100
# ./remesh ${path}water.obj ${path} water 500
# ./remesh ${path}tank.obj ${path} tank 100
# ---------------------------------------------------------------------------- #
# path='/Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/obj/Ren2015_no_float/'
# ./remesh ${path}wavemaker.obj ${path} wavemaker 100
# ./remesh ${path}water.obj ${path} water 500
# ./remesh ${path}tank.obj ${path} tank 100
# ---------------------------------------------------------------------------- #
# path='/Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/obj/tsukada2022_no_pool/'
# ./remesh ${path}water.obj ${path} water 1000
# ---------------------------------------------------------------------------- #
# MEMO: cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=remesh.cpp -DOUTPUT_NAME=remesh1
# path='/Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/obj/testALE/'
# ./remesh1 ${path}water.obj ${path} water_case1_ 1000 
# ./remesh2 ${path}water.obj ${path} water_case2_ 1000 
# ./remesh3 ${path}water.obj ${path} water_case3_ 1000 
# ./remesh4 ${path}water.obj ${path} water_case4_ 1000 
# ./remesh5 ${path}water.obj ${path} water_case5_ 1000 
# ./remesh5 ${path}cylinder.obj ${path} cylinder 100 &
# ./remesh5 ${path}cuboid.obj ${path} cuboid 100 &
# ./remesh5 ${path}tank.obj ${path} tank 100
# ---------------------------------------------------------------------------- #
path='/Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/obj/Chaplin2000/'
./remesh ${path}water.obj ${path} water 100 
# ---------------------------------------------------------------------------- #
# path='/Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/obj/Hadzic2005/'
# ./remesh ${path}water.obj ${path} water 1000 
# ./remesh ${path}wavemaker.obj ${path} wavemaker 500 
# ./remesh ${path}tank.obj ${path} tank 100 
# ./remesh ${path}float.obj ${path} float 100
# ---------------------------------------------------------------------------- #
# path='/Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/obj/Hadzic2005_24floats/'
# # from 1 to 24. from float1 to float24
# for i in {1..24}
# do
#     ./remesh ${path}float${i}.obj ${path} float${i}"_" 10
# done
# ./remesh ${path}water.obj ${path} water 1000
# ./remesh ${path}wavemaker.obj ${path} wavemaker 500 
# ./remesh ${path}tank.obj ${path} tank 100
