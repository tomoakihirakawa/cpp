#!/bin/sh

python3.11 input_generator.py -case Tanizawa1996 -mesh water_no_float0d09 -element pseudo_quad -wavemaker potential
python3.11 input_generator.py -case Tanizawa1996 -mesh water_no_float0d09 -element linear -wavemaker potential
python3.11 input_generator.py -case Tanizawa1996 -mesh water_no_float0d08 -element pseudo_quad -wavemaker potential
python3.11 input_generator.py -case Tanizawa1996 -mesh water_no_float0d08 -element linear -wavemaker potential
python3.11 input_generator.py -case Tanizawa1996 -mesh water_no_float0d07 -element pseudo_quad -wavemaker potential
python3.11 input_generator.py -case Tanizawa1996 -mesh water_no_float0d07 -element linear -wavemaker potential

./main ./input_files/Tanizawa1996_H0d05_L1d8_potential_water_no_float0d09_pseudo_quad
./main ./input_files/Tanizawa1996_H0d05_L1d8_potential_water_no_float0d08_pseudo_quad
./main ./input_files/Tanizawa1996_H0d05_L1d8_potential_water_no_float0d07_pseudo_quad

# ./main ./input_files/Tanizawa1996_H0d05_L1d8_potential_water_no_float0d09_linear
# ./main ./input_files/Tanizawa1996_H0d05_L1d8_potential_water_no_float0d08_linear
# ./main ./input_files/Tanizawa1996_H0d05_L1d8_potential_water_no_float0d07_linear

# python3.11 input_generator.py -case Tanizawa1996 -mesh water_no_float0d09 -element linear -wavemaker potential
# python3.11 input_generator.py -case Tanizawa1996 -mesh water_no_float0d08 -element linear -wavemaker potential
# python3.11 input_generator.py -case Tanizawa1996 -mesh water_no_float0d07 -element linear -wavemaker potential

# ./main ./input_files/Tanizawa1996_H0d05_L1d8_potential_water_no_float0d09_linear
# ./main ./input_files/Tanizawa1996_H0d05_L1d8_potential_water_no_float0d07_linear
# ./main ./input_files/Tanizawa1996_H0d05_L1d8_potential_water_no_float0d08_linear


# ./main ./input_files/Tanizawa1996_H0d05_L1d8gradPhiQuadElement__potential_water_no_float0d09
# ./main ./input_files/Tanizawa1996_H0d05_L1d8gradPhiQuadElement__potential_water_no_float0d08
# ./main ./input_files/Tanizawa1996_H0d05_L1d8gradPhiQuadElement__potential_water_no_float0d07

# ./main ./input_files/Tanizawa1996_H0d05_L1d8_quad_potential_water_no_float0d1
# ./main ./input_files/Tanizawa1996_H0d05_L1d8_quad_potential_water_no_float0d09
# ./main ./input_files/Tanizawa1996_H0d05_L1d8_quad_potential_water_no_float0d08
# ./main ./input_files/Tanizawa1996_H0d05_L1d8_quad_potential_water_no_float0d07
# ./main ./input_files/Tanizawa1996_H0d05_L1d8_quad_potential_water_no_float0d06

# ./main ./input_files/Tanizawa1996_H0d05_L1d8_potential_water_no_float0d1
# ./main ./input_files/Tanizawa1996_H0d05_L1d8_potential_water_no_float0d09
# ./main ./input_files/Tanizawa1996_H0d05_L1d8_potential_water_no_float0d08
# ./main ./input_files/Tanizawa1996_H0d05_L1d8_potential_water_no_float0d07
# ./main ./input_files/Tanizawa1996_H0d05_L1d8_potential_water_no_float0d06
# ./main ./input_files/Tanizawa1996_H0d05_L1d8_potential_water_no_float0d05

# ./main ./input_files/Tanizawa1996_H0d05_L1d8_potential_water0d1
# ./main ./input_files/Tanizawa1996_H0d05_L1d8_potential_water0d09
# ./main ./input_files/Tanizawa1996_H0d05_L1d8_potential_water0d08
# ./main ./input_files/Tanizawa1996_H0d05_L1d8_potential_water0d07
# ./main ./input_files/Tanizawa1996_H0d05_L1d8_potential_water0d06
# ./main ./input_files/Tanizawa1996_H0d05_L1d8_potential_water0d05

# ./main ./input_files/Tanizawa1996_H0d05_L1d8_potential
# ./main ./input_files/Tanizawa1996_H0d05_L2d7_potential

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