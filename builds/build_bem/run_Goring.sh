#!/bin/bash

python3.11 input_generator.py Goring1979 -m water0d1.obj -dt 0.05 -e pseudo_quad -o /Volumes/home/BEM/Goring1979/
python3.11 input_generator.py Goring1979 -m water0d1.obj -dt 0.05 -e linear -o /Volumes/home/BEM/Goring1979/
python3.11 input_generator.py Goring1979 -m water0d09.obj -dt 0.05 -e pseudo_quad -o /Volumes/home/BEM/Goring1979/
python3.11 input_generator.py Goring1979 -m water0d09.obj -dt 0.05 -e linear -o /Volumes/home/BEM/Goring1979/

# ./fast ./input_files/Goring1979_DT0d05_MESHwater0d09.obj_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
# ./fast ./input_files/Goring1979_DT0d05_MESHwater0d09.obj_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1
# ./fast ./input_files/Goring1979_DT0d05_MESHwater0d1.obj_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
# ./fast ./input_files/Goring1979_DT0d05_MESHwater0d1.obj_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1
