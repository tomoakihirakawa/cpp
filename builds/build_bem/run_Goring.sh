#!/bin/bash

python3.11 input_generator.py Goring1979 -dt 0.05 -s fast -e pseudo_quad -o /Volumes/home/BEM/Goring1979
python3.11 input_generator.py Goring1979 -dt 0.03 -s fast -e pseudo_quad -o /Volumes/home/BEM/Goring1979
python3.11 input_generator.py Goring1979 -dt 0.05 -s fast -o /Volumes/home/BEM/Goring1979
python3.11 input_generator.py Goring1979 -dt 0.03 -s fast -o /Volumes/home/BEM/Goring1979

./fast ./input_files/Goring1979Goring1979_DT0d03_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1_fast
./fast ./input_files/Goring1979Goring1979_DT0d03_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_fast
./fast ./input_files/Goring1979Goring1979_DT0d05_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1_fast
./fast ./input_files/Goring1979Goring1979_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_fast
