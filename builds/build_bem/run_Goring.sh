#!/bin/bash

python3.11 input_generator.py Goring1979 -m water0d09refined.obj -dt 0.03 -e pseudo_quad -ALE pseudo_quad -s 20250225
python3.11 input_generator.py Goring1979 -m water0d09refined.obj -dt 0.03 -e linear -ALE pseudo_quad -s 20250225
python3.11 input_generator.py Goring1979 -m water0d09refined.obj -dt 0.03 -e pseudo_quad -ALE linear -s 20250225
python3.11 input_generator.py Goring1979 -m water0d09refined.obj -dt 0.03 -e linear -ALE linear -s 20250225

python3.11 input_generator.py Goring1979 -m water0d06refined.obj -dt 0.03 -e pseudo_quad -ALE pseudo_quad -s 20250225
python3.11 input_generator.py Goring1979 -m water0d06refined.obj -dt 0.03 -e linear -ALE pseudo_quad -s 20250225
python3.11 input_generator.py Goring1979 -m water0d06refined.obj -dt 0.03 -e pseudo_quad -ALE linear -s 20250225
python3.11 input_generator.py Goring1979 -m water0d06refined.obj -dt 0.03 -e linear -ALE linear -s 20250225

# python3.11 input_generator.py Goring1979 -m water0d09refined.obj -dt 0.03 -e pseudo_quad -ALE pseudo_quad -o /Volumes/home/BEM/Goring1979/ -s 20250225
# python3.11 input_generator.py Goring1979 -m water0d09refined.obj -dt 0.03 -e linear -ALE pseudo_quad -o /Volumes/home/BEM/Goring1979/ -s 20250225
# python3.11 input_generator.py Goring1979 -m water0d08refined.obj -dt 0.03 -e pseudo_quad -ALE pseudo_quad -o /Volumes/home/BEM/Goring1979/ -s 20250225
# python3.11 input_generator.py Goring1979 -m water0d08refined.obj -dt 0.03 -e linear -ALE pseudo_quad -o /Volumes/home/BEM/Goring1979/ -s 20250225

# python3.11 input_generator.py Goring1979 -m water0d07refined.obj -dt 0.03 -e pseudo_quad -ALE linear -o /Volumes/home/BEM/Goring1979/ -s 20250225
# python3.11 input_generator.py Goring1979 -m water0d07refined.obj -dt 0.03 -e linear -ALE linear -o /Volumes/home/BEM/Goring1979/ -s 20250225
# python3.11 input_generator.py Goring1979 -m water0d07refined.obj -dt 0.03 -e pseudo_quad -ALE linear -o /Volumes/home/BEM/Goring1979/ -s 20250225
# python3.11 input_generator.py Goring1979 -m water0d07refined.obj -dt 0.03 -e linear -ALE linear -o /Volumes/home/BEM/Goring1979/ -s 20250225

./fast ./input_files/Goring1979_DT0d03_MESHwater0d06refined.obj_ELEMlinear_ALElinear_ALEPERIOD1_20250225
./fast ./input_files/Goring1979_DT0d03_MESHwater0d06refined.obj_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_20250225
./fast ./input_files/Goring1979_DT0d03_MESHwater0d06refined.obj_ELEMpseudo_quad_ALElinear_ALEPERIOD1_20250225
./fast ./input_files/Goring1979_DT0d03_MESHwater0d06refined.obj_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1_20250225

./fast ./input_files/Goring1979_DT0d03_MESHwater0d09refined.obj_ELEMlinear_ALElinear_ALEPERIOD1_20250225
./fast ./input_files/Goring1979_DT0d03_MESHwater0d09refined.obj_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_20250225
./fast ./input_files/Goring1979_DT0d03_MESHwater0d09refined.obj_ELEMpseudo_quad_ALElinear_ALEPERIOD1_20250225
./fast ./input_files/Goring1979_DT0d03_MESHwater0d09refined.obj_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1_20250225

# ./fast ./input_files/Goring1979_DT0d03_MESHwater0d06refined.obj_ELEMlinear_ALElinear_ALEPERIOD1_20250225
# ./fast ./input_files/Goring1979_DT0d03_MESHwater0d06refined.obj_ELEMpseudo_quad_ALElinear_ALEPERIOD1_20250225
# ./fast ./input_files/Goring1979_DT0d03_MESHwater0d06refined.obj_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_20250225
# ./fast ./input_files/Goring1979_DT0d03_MESHwater0d06refined.obj_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1_20250225

# ./fast ./input_files/Goring1979_DT0d03_MESHwater0d08refined.obj_ELEMlinear_ALElinear_ALEPERIOD1_20250225
# ./fast ./input_files/Goring1979_DT0d03_MESHwater0d08refined.obj_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_20250225
# ./fast ./input_files/Goring1979_DT0d03_MESHwater0d08refined.obj_ELEMpseudo_quad_ALElinear_ALEPERIOD1_20250225
# ./fast ./input_files/Goring1979_DT0d03_MESHwater0d08refined.obj_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1_20250225

# ./fast ./input_files/Goring1979_DT0d05_MESHwater0d08.obj_ELEMlinear_ALElinear_ALEPERIOD1_20250225
# ./fast ./input_files/Goring1979_DT0d05_MESHwater0d08.obj_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_20250225
# ./fast ./input_files/Goring1979_DT0d05_MESHwater0d08.obj_ELEMpseudo_quad_ALElinear_ALEPERIOD1_20250225
# ./fast ./input_files/Goring1979_DT0d05_MESHwater0d08.obj_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1_20250225

# ./fast ./input_files/Goring1979_DT0d05_MESHwater0d07.obj_ELEMlinear_ALElinear_ALEPERIOD1_20250225
# ./fast ./input_files/Goring1979_DT0d05_MESHwater0d07.obj_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_20250225
# ./fast ./input_files/Goring1979_DT0d05_MESHwater0d07.obj_ELEMpseudo_quad_ALElinear_ALEPERIOD1_20250225
# ./fast ./input_files/Goring1979_DT0d05_MESHwater0d07.obj_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1_20250225
