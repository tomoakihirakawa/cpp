#!/bin/bash

outputdir=/Volumes/home/BEM/benchmark202411
case=Goring1979
mesh_array=(water100 water150 water200)
dt_array=(0.02 0.04 0.06)

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

./main ./input_files/Goring1979_DT0d06_ELEMlinear_ALElinear_ALEPERIOD1
./main ./input_files/Goring1979_DT0d06_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
./main ./input_files/Goring1979_DT0d06_ELEMpseudo_quad_ALElinear_ALEPERIOD1
./main ./input_files/Goring1979_DT0d06_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1

./main ./input_files/Goring1979_DT0d04_ELEMlinear_ALElinear_ALEPERIOD1
./main ./input_files/Goring1979_DT0d04_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
./main ./input_files/Goring1979_DT0d04_ELEMpseudo_quad_ALElinear_ALEPERIOD1
./main ./input_files/Goring1979_DT0d04_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1

./main ./input_files/Goring1979_DT0d02_ELEMlinear_ALElinear_ALEPERIOD1
./main ./input_files/Goring1979_DT0d02_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
./main ./input_files/Goring1979_DT0d02_ELEMpseudo_quad_ALElinear_ALEPERIOD1
./main ./input_files/Goring1979_DT0d02_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1
