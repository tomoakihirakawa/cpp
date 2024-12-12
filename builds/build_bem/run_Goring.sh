#!/bin/bash

outputdir=/Volumes/home/BEM/benchmark202411
case=Goring1979
mesh_array=(water0d1 water0d09 water0d08)
dt_array=(0.05 0.1)

wavemaker=Goring1979
element_array=(pseudo_quad linear)
ALE_array=(pseudo_quad linear)

for ALE in ${ALE_array[@]}; do
    for element in ${element_array[@]}; do
        for dt in ${dt_array[@]}; do
            for mesh in ${mesh_array[@]}; do
                python3.11 input_generator.py -case ${case} -mesh ${mesh} -element ${element} -wavemaker ${wavemaker} -dt ${dt} -ALE ${ALE} -outputdir ${outputdir} -ALEPERIOD 1
            done
        done
    done
done

# ./main ./input_files/Goring1979_MESHwater0d1_DT0d1_ELEMlinear_ALElinear_ALEPERIOD1
# ./main ./input_files/Goring1979_MESHwater0d1_DT0d1_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
# ./main ./input_files/Goring1979_MESHwater0d1_DT0d1_ELEMpseudo_quad_ALElinear_ALEPERIOD1
# ./main ./input_files/Goring1979_MESHwater0d1_DT0d1_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1

# ./main ./input_files/Goring1979_MESHwater0d08_DT0d1_ELEMlinear_ALElinear_ALEPERIOD1
# ./main ./input_files/Goring1979_MESHwater0d08_DT0d1_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
# ./main ./input_files/Goring1979_MESHwater0d08_DT0d1_ELEMpseudo_quad_ALElinear_ALEPERIOD1
# ./main ./input_files/Goring1979_MESHwater0d08_DT0d1_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1

# ./main ./input_files/Goring1979_MESHwater0d09_DT0d1_ELEMlinear_ALElinear_ALEPERIOD1
# ./main ./input_files/Goring1979_MESHwater0d09_DT0d1_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
# ./main ./input_files/Goring1979_MESHwater0d09_DT0d1_ELEMpseudo_quad_ALElinear_ALEPERIOD1
# ./main ./input_files/Goring1979_MESHwater0d09_DT0d1_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1

# ./main ./input_files/Goring1979_MESHwater0d1_DT0d05_ELEMlinear_ALElinear_ALEPERIOD1
# ./main ./input_files/Goring1979_MESHwater0d1_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
# ./main ./input_files/Goring1979_MESHwater0d1_DT0d05_ELEMpseudo_quad_ALElinear_ALEPERIOD1
# ./main ./input_files/Goring1979_MESHwater0d1_DT0d05_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1

# ./main ./input_files/Goring1979_MESHwater0d08_DT0d05_ELEMlinear_ALElinear_ALEPERIOD1
# ./main ./input_files/Goring1979_MESHwater0d08_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
# ./main ./input_files/Goring1979_MESHwater0d08_DT0d05_ELEMpseudo_quad_ALElinear_ALEPERIOD1
# ./main ./input_files/Goring1979_MESHwater0d08_DT0d05_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1

# ./main ./input_files/Goring1979_MESHwater0d09_DT0d05_ELEMlinear_ALElinear_ALEPERIOD1
# ./main ./input_files/Goring1979_MESHwater0d09_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
# ./main ./input_files/Goring1979_MESHwater0d09_DT0d05_ELEMpseudo_quad_ALElinear_ALEPERIOD1
# ./main ./input_files/Goring1979_MESHwater0d09_DT0d05_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1
