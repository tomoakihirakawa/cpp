#!/bin/sh

# outputdir=/Volumes/home/BEM/benchmark202404

# set home directory
outputdir=/Users/tomoaki

case=Fredriksen2015

# step 0.05 from 0.6 to 1.0

for T in 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0; do
    python3.11 input_generator.py ${case} -dt 0.05 -T ${T} -o ${outputdir}/BEM/benchmark20250117_Fredriksen2015 -e pseudo_quad
done

# ./main ./input_files/Fredriksen2015_WS0d01667_T0d9_DT0d05_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1
# ./main ./input_files/Fredriksen2015_WS0d01667_T0d95_DT0d05_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1
# ./main ./input_files/Fredriksen2015_WS0d01667_T1d0_DT0d05_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1
# ./main ./input_files/Fredriksen2015_WS0d01667_T0d6_DT0d05_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1
# ./main ./input_files/Fredriksen2015_WS0d01667_T0d65_DT0d05_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1
# ./main ./input_files/Fredriksen2015_WS0d01667_T0d7_DT0d05_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1
# ./main ./input_files/Fredriksen2015_WS0d01667_T0d75_DT0d05_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1
# ./main ./input_files/Fredriksen2015_WS0d01667_T0d8_DT0d05_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1
# ./main ./input_files/Fredriksen2015_WS0d01667_T0d85_DT0d05_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1
