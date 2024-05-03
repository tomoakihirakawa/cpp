#!/bin/sh

# outputdir=/Volumes/home/BEM/benchmark202404
outputdir=${HOME}/BEM

case=Tanizawa1996
mesh_array=("water_no_float0d08" "water_no_float0d07")
dt_array=(0.05 0.1)
H_array=(0.04 0.08)

# outputdir=/Volumes/home/BEM/benchmark202404Palm2016
# case=Palm2016
# mesh_array=(water_mod)
# dt_array=(0.05 0.1)
# H_array=(0.04 0.2)

wavemaker=potential
element_array=(pseudo_quad linear)
ALE_array=(pseudo_quad linear)

for ALE in ${ALE_array[@]};do
    for element in ${element_array[@]};do
        for dt in ${dt_array[@]};do
            for H in ${H_array[@]};do
                for mesh in ${mesh_array[@]};do
                    python3.11 input_generator.py -case ${case} -mesh ${mesh} -element ${element} -wavemaker ${wavemaker} -dt ${dt} -H ${H} -ALE ${ALE} -outputdir ${outputdir}
                done
            done
        done
    done
done

# ./main ./input_files/Palm2016_H0d2_T1d0_MESHwater_mod_DT0d1_ELEMpseudo_quad_ALEpseudo_quad

# ./main ./input_files/${case}_H0d3_L1d8_WAVEflap_MESHwater_no_float0d08_DT0d05_ELEMpseudo_quad_ALEpseudo_quad

# ----------------------------------- H005 ----------------------------------- #

# mac studio
./main ./input_files/${case}_H0d05_L1d8_WAVEflap_MESHwater_no_float0d08_DT0d05_ELEMpseudo_quad_ALEpseudo_quad
./main ./input_files/${case}_H0d05_L1d8_WAVEflap_MESHwater_no_float0d08_DT0d05_ELEMlinear_ALEpseudo_quad
./main ./input_files/${case}_H0d05_L1d8_WAVEflap_MESHwater_no_float0d08_DT0d05_ELEMpseudo_quad_ALElinear
./main ./input_files/${case}_H0d05_L1d8_WAVEflap_MESHwater_no_float0d08_DT0d05_ELEMlinear_ALElinear

# ./main ./input_files/${case}_H0d05_L1d8_WAVEflap_MESHwater_no_float0d08_DT0d03_ELEMpseudo_quad_ALEpseudo_quad # prandtl 142
# ./main ./input_files/${case}_H0d05_L1d8_WAVEflap_MESHwater_no_float0d08_DT0d03_ELEMlinear_ALEpseudo_quad # Daiki 140
# ./main ./input_files/${case}_H0d05_L1d8_WAVEflap_MESHwater_no_float0d08_DT0d03_ELEMpseudo_quad_ALElinear #macbook 
# ./main ./input_files/${case}_H0d05_L1d8_WAVEflap_MESHwater_no_float0d08_DT0d03_ELEMlinear_ALElinear #ishiwaka 141

# ------------------------------------ H01 ----------------------------------- #

# ./main ./input_files/${case}_H0d1_L1d8_WAVEflap_MESHwater_no_float0d08_DT0d05_ELEMpseudo_quad_ALEpseudo_quad
# ./main ./input_files/${case}_H0d1_L1d8_WAVEflap_MESHwater_no_float0d08_DT0d05_ELEMlinear_ALEpseudo_quad
# ./main ./input_files/${case}_H0d1_L1d8_WAVEflap_MESHwater_no_float0d08_DT0d05_ELEMpseudo_quad_ALElinear
# ./main ./input_files/${case}_H0d1_L1d8_WAVEflap_MESHwater_no_float0d08_DT0d05_ELEMlinear_ALElinear

# ./main ./input_files/${case}_H0d1_L1d8_WAVEflap_MESHwater_no_float0d08_DT0d03_ELEMpseudo_quad_ALEpseudo_quad # prandtl 142
# ./main ./input_files/${case}_H0d1_L1d8_WAVEflap_MESHwater_no_float0d08_DT0d03_ELEMlinear_ALEpseudo_quad# Daiki 140
# ./main ./input_files/${case}_H0d1_L1d8_WAVEflap_MESHwater_no_float0d08_DT0d03_ELEMpseudo_quad_ALElinear #macbook 
# ./main ./input_files/${case}_H0d1_L1d8_WAVEflap_MESHwater_no_float0d08_DT0d03_ELEMlinear_ALElinear #ishiwaka 141
