#!/bin/sh

# for T in 7d5 5d0 5d5 6d0 6d5 7d0 8d0 8d5 9d0 9d5 10d0; do
#   file="./input_files/moon_pool_large_a0d8_T${T}_h80_modified_mesh"
#   ./main ${file}
#   rsync -v ${file} Kelvin@10.0.1.14:~/BEM/
# done

./main ./input_files/Ren2015T1d2_a0d06
./main ./input_files/Ren2015T1d2_a0d04