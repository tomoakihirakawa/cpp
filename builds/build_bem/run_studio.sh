#!/bin/sh

outputdir=/Volumes/home/BEM/benchmark202405
# outputdir=/Volumes/home/BEM/benchmark202404Palm2016

# outputdir=${HOME}/BEM
case=Tanizawa1996
mesh_array=(water_no_float0d08 water_no_float0d07 water_no_float0d06)
dt_array=(0.05 0.025)
H_array=(0.05 0.1)
ALEPERIOD_array=(1 3)

# case=Palm2016
# mesh_array=(water_mod)
# dt_array=(0.05 0.1)
# H_array=(0.04 0.2)

wavemaker=potential
element_array=(pseudo_quad linear)
ALE_array=(pseudo_quad linear)

for ALEPERIOD in ${ALEPERIOD_array[@]};do
    for ALE in ${ALE_array[@]};do
        for element in ${element_array[@]};do
            for dt in ${dt_array[@]};do
                for H in ${H_array[@]};do
                    for mesh in ${mesh_array[@]};do
                        python3.11 input_generator.py -case ${case} -mesh ${mesh} -element ${element} -wavemaker ${wavemaker} -dt ${dt} -H ${H} -ALE ${ALE} -outputdir ${outputdir} -ALEPERIOD ${ALEPERIOD}
                    done
                done
            done
        done
    done
done

./main ./input_files/Tanizawa1996_H0d05_L1d8_WAVEpotential_MESHwater_no_float0d08_DT0d05_ELEMlinear_ALElinear_ALEPERIOD1
./main ./input_files/Tanizawa1996_H0d05_L1d8_WAVEpotential_MESHwater_no_float0d07_DT0d05_ELEMlinear_ALElinear_ALEPERIOD1
./main ./input_files/Tanizawa1996_H0d05_L1d8_WAVEpotential_MESHwater_no_float0d06_DT0d05_ELEMlinear_ALElinear_ALEPERIOD1
# 
./main ./input_files/Tanizawa1996_H0d1_L1d8_WAVEpotential_MESHwater_no_float0d08_DT0d05_ELEMlinear_ALElinear_ALEPERIOD1
./main ./input_files/Tanizawa1996_H0d1_L1d8_WAVEpotential_MESHwater_no_float0d07_DT0d05_ELEMlinear_ALElinear_ALEPERIOD1
./main ./input_files/Tanizawa1996_H0d1_L1d8_WAVEpotential_MESHwater_no_float0d06_DT0d05_ELEMlinear_ALElinear_ALEPERIOD1
./main ./input_files/Tanizawa1996_H0d1_L1d8_WAVEpotential_MESHwater_no_float0d08_DT0d05_ELEMlinear_ALElinear_ALEPERIOD3
./main ./input_files/Tanizawa1996_H0d1_L1d8_WAVEpotential_MESHwater_no_float0d07_DT0d05_ELEMlinear_ALElinear_ALEPERIOD3
./main ./input_files/Tanizawa1996_H0d1_L1d8_WAVEpotential_MESHwater_no_float0d06_DT0d05_ELEMlinear_ALElinear_ALEPERIOD3
# 
./main ./input_files/Tanizawa1996_H0d1_L1d8_WAVEpotential_MESHwater_no_float0d08_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
./main ./input_files/Tanizawa1996_H0d1_L1d8_WAVEpotential_MESHwater_no_float0d07_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
./main ./input_files/Tanizawa1996_H0d1_L1d8_WAVEpotential_MESHwater_no_float0d06_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
./main ./input_files/Tanizawa1996_H0d1_L1d8_WAVEpotential_MESHwater_no_float0d08_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD3
./main ./input_files/Tanizawa1996_H0d1_L1d8_WAVEpotential_MESHwater_no_float0d07_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD3
./main ./input_files/Tanizawa1996_H0d1_L1d8_WAVEpotential_MESHwater_no_float0d06_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD3
# 
./main ./input_files/Tanizawa1996_H0d1_L1d8_WAVEpotential_MESHwater_no_float0d08_DT0d05_ELEMpseudo_quad_ALElinear_ALEPERIOD1
./main ./input_files/Tanizawa1996_H0d1_L1d8_WAVEpotential_MESHwater_no_float0d07_DT0d05_ELEMpseudo_quad_ALElinear_ALEPERIOD1
./main ./input_files/Tanizawa1996_H0d1_L1d8_WAVEpotential_MESHwater_no_float0d06_DT0d05_ELEMpseudo_quad_ALElinear_ALEPERIOD1
./main ./input_files/Tanizawa1996_H0d1_L1d8_WAVEpotential_MESHwater_no_float0d08_DT0d05_ELEMpseudo_quad_ALElinear_ALEPERIOD3
./main ./input_files/Tanizawa1996_H0d1_L1d8_WAVEpotential_MESHwater_no_float0d07_DT0d05_ELEMpseudo_quad_ALElinear_ALEPERIOD3
./main ./input_files/Tanizawa1996_H0d1_L1d8_WAVEpotential_MESHwater_no_float0d06_DT0d05_ELEMpseudo_quad_ALElinear_ALEPERIOD3

# macstudio
# ./main ./input_files/Tanizawa1996_H0d1_L1d8_WAVEpotential_MESHwater_no_float0d08_DT0d05_ELEMlinear_ALElinear_ALEPERIOD3
# ./main ./input_files/Tanizawa1996_H0d1_L1d8_WAVEpotential_MESHwater_no_float0d08_DT0d05_ELEMlinear_ALElinear_ALEPERIOD1
# ./main ./input_files/Tanizawa1996_H0d1_L1d8_WAVEpotential_MESHwater_no_float0d08_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD3
# ./main ./input_files/Tanizawa1996_H0d1_L1d8_WAVEpotential_MESHwater_no_float0d08_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1

# ./main ./input_files/Tanizawa1996_H0d1_L1d8_WAVEpotential_MESHwater_no_float0d08_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD3
# ./main ./input_files/Tanizawa1996_H0d1_L1d8_WAVEpotential_MESHwater_no_float0d07_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD3
# ./main ./input_files/Tanizawa1996_H0d1_L1d8_WAVEpotential_MESHwater_no_float0d06_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD3

# ./main ./input_files/Tanizawa1996_H0d1_L1d8_WAVEpotential_MESHwater_no_float0d08_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
# ./main ./input_files/Tanizawa1996_H0d1_L1d8_WAVEpotential_MESHwater_no_float0d07_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
# ./main ./input_files/Tanizawa1996_H0d1_L1d8_WAVEpotential_MESHwater_no_float0d06_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1

# # 以下はできれば

# ./main ./input_files/Tanizawa1996_H0d1_L1d8_WAVEpotential_MESHwater_no_float0d07_DT0d05_ELEMlinear_ALElinear_ALEPERIOD3
# ./main ./input_files/Tanizawa1996_H0d1_L1d8_WAVEpotential_MESHwater_no_float0d07_DT0d05_ELEMlinear_ALElinear_ALEPERIOD1

# ./main ./input_files/Tanizawa1996_H0d1_L1d8_WAVEpotential_MESHwater_no_float0d07_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD3
# ./main ./input_files/Tanizawa1996_H0d1_L1d8_WAVEpotential_MESHwater_no_float0d07_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
