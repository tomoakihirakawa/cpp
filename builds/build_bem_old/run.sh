#!/bin/sh

./main ./input_files/Palm2016_with_mooring
./main ./input_files/Palm2016_without_mooring
# ./main ./input_files/Kramer2021_H00d09_small_meshA_ALE3
# ./main ./input_files/Kramer2021_H00d09_small_meshB_ALE3
# ./main ./input_files/Kramer2021_H00d09_small_meshC_ALE3
# ./main ./input_files/Kramer2021_H00d09_small_meshF_ALE1
./main ./input_files/Kramer2021_H00d09_small_meshE_ALE1
./main ./input_files/Kramer2021_H00d09_small_meshD_ALE1

# ./main ./input_files/Cheng2018meshB_H0d04_T1d2_piston_positive
# ./main ./input_files/Cheng2018meshF_H0d04_T1d2_piston_positive

# ./main ./input_files/Cheng2018meshC_H0d1_T1d2_piston_positive
# ./main ./input_files/Cheng2018meshD_H0d1_T1d2_piston_positive
# ./main ./input_files/Cheng2018meshE_H0d1_T1d2_piston_positive
# ./main ./input_files/Cheng2018meshF_H0d1_T1d2_piston_positive
# ./main ./input_files/Cheng2018meshA_H0d1_T1d2_piston_positive
# ./main ./input_files/Cheng2018meshB_H0d1_T1d2_piston_positive

# --------------------------------- moon pool -------------------------------- #

# for barge in large no; do
# for T in 7d5 5d0 5d5 6d0 6d5 7d0 8d0 8d5 9d0 9d5 10d0; do
#   file="./input_files/moon_pool_${barge}_a0d8_T${T}_h80"
#   ./main ${file}
# #   rsync -v ${file} Kelvin@10.0.1.14:~/BEM/
# done
# done

# -------------------------------- Kramer2021 -------------------------------- #

# ./main ./input_files/Kramer2021_H00d03_small

# -------------------------------- Hadzic2005 -------------------------------- #

# ./main ./input_files/Hadzic2005