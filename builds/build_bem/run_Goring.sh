#!/bin/bash

outputdir=/Volumes/home/BEM/benchmark202411
case=Goring1979
mesh_array=(water120 water140 water160)
dt_array=(0.03 0.06)

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

# ./main ./input_files/Goring1979_MESHwater120_DT0d06_ELEMlinear_ALElinear_ALEPERIOD1
# ./main ./input_files/Goring1979_MESHwater120_DT0d06_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
# ./main ./input_files/Goring1979_MESHwater120_DT0d06_ELEMpseudo_quad_ALElinear_ALEPERIOD1
./main ./input_files/Goring1979_MESHwater120_DT0d06_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1

# ./main ./input_files/Goring1979_MESHwater140_DT0d06_ELEMlinear_ALElinear_ALEPERIOD1
# ./main ./input_files/Goring1979_MESHwater140_DT0d06_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
# ./main ./input_files/Goring1979_MESHwater140_DT0d06_ELEMpseudo_quad_ALElinear_ALEPERIOD1
./main ./input_files/Goring1979_MESHwater140_DT0d06_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1

# ./main ./input_files/Goring1979_MESHwater160_DT0d06_ELEMlinear_ALElinear_ALEPERIOD1
# ./main ./input_files/Goring1979_MESHwater160_DT0d06_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
# ./main ./input_files/Goring1979_MESHwater160_DT0d06_ELEMpseudo_quad_ALElinear_ALEPERIOD1
./main ./input_files/Goring1979_MESHwater160_DT0d06_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1

# ./main ./input_files/Goring1979_MESHwater120_DT0d03_ELEMlinear_ALElinear_ALEPERIOD1
# ./main ./input_files/Goring1979_MESHwater120_DT0d03_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
# ./main ./input_files/Goring1979_MESHwater120_DT0d03_ELEMpseudo_quad_ALElinear_ALEPERIOD1
./main ./input_files/Goring1979_MESHwater120_DT0d03_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1

# ./main ./input_files/Goring1979_MESHwater140_DT0d03_ELEMlinear_ALElinear_ALEPERIOD1
# ./main ./input_files/Goring1979_MESHwater140_DT0d03_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
# ./main ./input_files/Goring1979_MESHwater140_DT0d03_ELEMpseudo_quad_ALElinear_ALEPERIOD1
./main ./input_files/Goring1979_MESHwater140_DT0d03_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1

# ./main ./input_files/Goring1979_MESHwater160_DT0d03_ELEMlinear_ALElinear_ALEPERIOD1
# ./main ./input_files/Goring1979_MESHwater160_DT0d03_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
# ./main ./input_files/Goring1979_MESHwater160_DT0d03_ELEMpseudo_quad_ALElinear_ALEPERIOD1
./main ./input_files/Goring1979_MESHwater160_DT0d03_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1
