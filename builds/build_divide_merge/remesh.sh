#! /bin/bash

# ./remesh ~/Dropbox/code/cpp/obj/Ren2015/wavemaker.obj ~/Dropbox/code/cpp/obj/Ren2015 wavemaker 100
# ./remesh ~/Dropbox/code/cpp/obj/Ren2015/float.obj ~/Dropbox/code/cpp/obj/Ren2015 float 100
# ./remesh ~/Dropbox/code/cpp/obj/Ren2015/water.obj ~/Dropbox/code/cpp/obj/Ren2015 water 400
# ./remesh ~/Dropbox/code/cpp/obj/Ren2015/water.obj ~/Dropbox/code/cpp/obj/Ren2015 water 100
# ./remesh ~/Dropbox/code/cpp/obj/Ren2015/tank.obj ~/Dropbox/code/cpp/obj/Ren2015 tank 100

# define string variable

path='/Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/obj/Ren2015_multiple/'

# ./remesh ${path}float_a.obj ${path} float_a 100
# ./remesh ${path}float_b.obj ${path} float_b 100
# ./remesh ${path}float_c.obj ${path} float_c 100
# ./remesh ${path}float_d.obj ${path} float_d 100
# ./remesh ${path}float_e.obj ${path} float_e 100
# ./remesh ${path}float_f.obj ${path} float_f 100
# ./remesh ${path}float_g.obj ${path} float_g 100

# ./remesh ${path}wavemaker.obj ${path} wavemaker 100
./remesh ${path}water.obj ${path} water 500
./remesh ${path}tank.obj ${path} tank 100