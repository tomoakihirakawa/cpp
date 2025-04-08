#!/bin/sh

# これは海洋開発論文集のためのベンチマークです．

cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_NAME=main.cpp -DOUTPUT_NAME=fast

python3.11 input_generator.py Kramer2021 -dt 0.01 -s H00d03 -o /Volumes/home/BEM/benchmark20250404
./fast ./input_files/Kramer2021_DT0d01_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_H00d03

python3.11 input_generator.py Tanizawa1996 -dt 0.1 -H 0.05 -L 1.8 -e linear -w potential -o /Volumes/home/BEM/benchmark20250404
python3.11 input_generator.py Tanizawa1996 -dt 0.05 -H 0.05 -L 1.8 -e linear -w potential -o /Volumes/home/BEM/benchmark20250404
./fast ./input_files/Tanizawa1996_H0d05_L1d8_DT0d1_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
./fast ./input_files/Tanizawa1996_H0d05_L1d8_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1

python3.11 input_generator.py DeepCWind -dt 0.5 -H 2 -o /Volumes/home/BEM/benchmark20250401 -s 4floats
./fast ./input_files/DeepCWind_DT0d5_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_4floats
