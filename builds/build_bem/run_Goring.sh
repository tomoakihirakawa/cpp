#!/bin/bash

python3.11 input_generator.py Goring1979 -m water0d09.obj -dt 0.05 -e pseudo_quad -ALE pseudo_quad -o /Volumes/home/BEM/Goring1979/
python3.11 input_generator.py Goring1979 -m water0d09.obj -dt 0.05 -e linear -ALE pseudo_quad -o /Volumes/home/BEM/Goring1979/
python3.11 input_generator.py Goring1979 -m water0d08.obj -dt 0.05 -e pseudo_quad -ALE pseudo_quad -o /Volumes/home/BEM/Goring1979/
python3.11 input_generator.py Goring1979 -m water0d08.obj -dt 0.05 -e linear -ALE pseudo_quad -o /Volumes/home/BEM/Goring1979/
python3.11 input_generator.py Goring1979 -m water0d07.obj -dt 0.05 -e pseudo_quad -ALE pseudo_quad -o /Volumes/home/BEM/Goring1979/
python3.11 input_generator.py Goring1979 -m water0d07.obj -dt 0.05 -e linear -ALE pseudo_quad -o /Volumes/home/BEM/Goring1979/

python3.11 input_generator.py Goring1979 -m water0d09.obj -dt 0.05 -e pseudo_quad -ALE linear -o /Volumes/home/BEM/Goring1979/
python3.11 input_generator.py Goring1979 -m water0d09.obj -dt 0.05 -e linear -ALE linear -o /Volumes/home/BEM/Goring1979/
python3.11 input_generator.py Goring1979 -m water0d08.obj -dt 0.05 -e pseudo_quad -ALE linear -o /Volumes/home/BEM/Goring1979/
python3.11 input_generator.py Goring1979 -m water0d08.obj -dt 0.05 -e linear -ALE linear -o /Volumes/home/BEM/Goring1979/
python3.11 input_generator.py Goring1979 -m water0d07.obj -dt 0.05 -e pseudo_quad -ALE linear -o /Volumes/home/BEM/Goring1979/
python3.11 input_generator.py Goring1979 -m water0d07.obj -dt 0.05 -e linear -ALE linear -o /Volumes/home/BEM/Goring1979/

# ./fast ./input_files/Goring1979_DT0d05_MESHwater0d09.obj_ELEMlinear_ALElinear_ALEPERIOD1
# ./fast ./input_files/Goring1979_DT0d05_MESHwater0d09.obj_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
# ./fast ./input_files/Goring1979_DT0d05_MESHwater0d09.obj_ELEMpseudo_quad_ALElinear_ALEPERIOD1
# ./fast ./input_files/Goring1979_DT0d05_MESHwater0d09.obj_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1

# ./fast ./input_files/Goring1979_DT0d05_MESHwater0d08.obj_ELEMlinear_ALElinear_ALEPERIOD1
# ./fast ./input_files/Goring1979_DT0d05_MESHwater0d08.obj_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
# ./fast ./input_files/Goring1979_DT0d05_MESHwater0d08.obj_ELEMpseudo_quad_ALElinear_ALEPERIOD1
# ./fast ./input_files/Goring1979_DT0d05_MESHwater0d08.obj_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1

# ./fast ./input_files/Goring1979_DT0d05_MESHwater0d07.obj_ELEMlinear_ALElinear_ALEPERIOD1
# ./fast ./input_files/Goring1979_DT0d05_MESHwater0d07.obj_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
# ./fast ./input_files/Goring1979_DT0d05_MESHwater0d07.obj_ELEMpseudo_quad_ALElinear_ALEPERIOD1
# ./fast ./input_files/Goring1979_DT0d05_MESHwater0d07.obj_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1
