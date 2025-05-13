#!/bin/sh

# これは海洋開発論文集のためのベンチマークです．

# cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_NAME=main.cpp -DOUTPUT_NAME=fast

# python3.11 input_generator.py Kramer2021 -dt 0.01 -s H00d03 -o /Volumes/home/BEM/benchmark20250404
# ./fast ./input_files/Kramer2021_DT0d01_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_H00d03

# python3.11 input_generator.py Tanizawa1996 -dt 0.1 -H 0.05 -L 1.8 -e linear -w potential -o /Volumes/home/BEM/benchmark20250408
# python3.11 input_generator.py Tanizawa1996 -dt 0.05 -H 0.05 -L 1.8 -e linear -w potential -o /Volumes/home/BEM/benchmark20250404
# ./fast ./input_files/Tanizawa1996_H0d05_L1d8_DT0d1_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
# ./fast ./input_files/Tanizawa1996_H0d05_L1d8_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1

# python3.11 input_generator.py DeepCWind -dt 0.5 -H 2 -o /Volumes/home/BEM/benchmark20250401 -s 4floats
# ./fast ./input_files/DeepCWind_DT0d5_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_4floats

python3 input_generator.py Ruehl2016 -m water0d15.obj -dt 0.1 -e pseudo_quad -ALE pseudo_quad -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501 -s Trial04
python3 input_generator.py Ruehl2016 -m water0d15.obj -dt 0.1 -e pseudo_quad -ALE linear -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501 -s Trial04
python3 input_generator.py Ruehl2016 -m water0d15.obj -dt 0.1 -e linear -ALE pseudo_quad -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501 -s Trial04
python3 input_generator.py Ruehl2016 -m water0d15.obj -dt 0.1 -e linear -ALE linear -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501 -s Trial04

python3 input_generator.py Ruehl2016 -m water0d15.obj -dt 0.1 -e pseudo_quad -ALE pseudo_quad -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501 -s Trial05
python3 input_generator.py Ruehl2016 -m water0d15.obj -dt 0.1 -e pseudo_quad -ALE linear -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501 -s Trial05
python3 input_generator.py Ruehl2016 -m water0d15.obj -dt 0.1 -e linear -ALE pseudo_quad -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501 -s Trial05
python3 input_generator.py Ruehl2016 -m water0d15.obj -dt 0.1 -e linear -ALE linear -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501 -s Trial05

python3 input_generator.py Ruehl2016 -m water0d15.obj -dt 0.1 -e pseudo_quad -ALE pseudo_quad -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501 -s Trial06
python3 input_generator.py Ruehl2016 -m water0d15.obj -dt 0.1 -e pseudo_quad -ALE linear -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501 -s Trial06
python3 input_generator.py Ruehl2016 -m water0d15.obj -dt 0.1 -e linear -ALE pseudo_quad -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501 -s Trial06
python3 input_generator.py Ruehl2016 -m water0d15.obj -dt 0.1 -e linear -ALE linear -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501 -s Trial06

python3 input_generator.py Ruehl2016 -m water0d15.obj -dt 0.1 -e pseudo_quad -ALE pseudo_quad -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501 -s Trial09
python3 input_generator.py Ruehl2016 -m water0d15.obj -dt 0.1 -e pseudo_quad -ALE linear -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501 -s Trial09
python3 input_generator.py Ruehl2016 -m water0d15.obj -dt 0.1 -e linear -ALE pseudo_quad -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501 -s Trial09
python3 input_generator.py Ruehl2016 -m water0d15.obj -dt 0.1 -e linear -ALE linear -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501 -s Trial09

python3 input_generator.py Ruehl2016 -m water0d15.obj -dt 0.1 -e pseudo_quad -ALE pseudo_quad -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501 -s Trial08
python3 input_generator.py Ruehl2016 -m water0d15.obj -dt 0.1 -e pseudo_quad -ALE linear -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501 -s Trial08
python3 input_generator.py Ruehl2016 -m water0d15.obj -dt 0.1 -e linear -ALE pseudo_quad -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501 -s Trial08
python3 input_generator.py Ruehl2016 -m water0d15.obj -dt 0.1 -e linear -ALE linear -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501 -s Trial08

python3 input_generator.py Ruehl2016 -m water0d15.obj -dt 0.1 -e pseudo_quad -ALE pseudo_quad -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501 -s Trial07
python3 input_generator.py Ruehl2016 -m water0d15.obj -dt 0.1 -e pseudo_quad -ALE linear -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501 -s Trial07
python3 input_generator.py Ruehl2016 -m water0d15.obj -dt 0.1 -e linear -ALE pseudo_quad -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501 -s Trial07
python3 input_generator.py Ruehl2016 -m water0d15.obj -dt 0.1 -e linear -ALE linear -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501 -s Trial07

./fast ./input_files/Ruehl2016_DT0d1_MESHwater0d15.obj_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1_Trial09
./fast ./input_files/Ruehl2016_DT0d1_MESHwater0d15.obj_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_Trial09
./fast ./input_files/Ruehl2016_DT0d1_MESHwater0d15.obj_ELEMlinear_ALElinear_ALEPERIOD1_Trial09
./fast ./input_files/Ruehl2016_DT0d1_MESHwater0d15.obj_ELEMpseudo_quad_ALElinear_ALEPERIOD 1_Trial09

# 1. バケツがfaceを保持していないか，この場合少しのジャンプだろう（おそらくこっちか）
# 2. 保持していたとしても，バケツ自体を取り出せていないかもしれない．この場合大きなジャンプが生じていただろう
# 結果で確認．

./fast ./input_files/Ruehl2016_DT0d1_MESHwater0d15.obj_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1_Trial08
./fast ./input_files/Ruehl2016_DT0d1_MESHwater0d15.obj_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_Trial08
./fast ./input_files/Ruehl2016_DT0d1_MESHwater0d15.obj_ELEMlinear_ALElinear_ALEPERIOD1_Trial08
./fast ./input_files/Ruehl2016_DT0d1_MESHwater0d15.obj_ELEMpseudo_quad_ALElinear_ALEPERIOD1_Trial08

./fast ./input_files/Ruehl2016_DT0d1_MESHwater0d15.obj_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1_Trial07
./fast ./input_files/Ruehl2016_DT0d1_MESHwater0d15.obj_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_Trial07
./fast ./input_files/Ruehl2016_DT0d1_MESHwater0d15.obj_ELEMlinear_ALElinear_ALEPERIOD1_Trial07
./fast ./input_files/Ruehl2016_DT0d1_MESHwater0d15.obj_ELEMpseudo_quad_ALElinear_ALEPERIOD1_Trial07

./fast ./input_files/Ruehl2016_DT0d1_MESHwater0d15.obj_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1_Trial06
./fast ./input_files/Ruehl2016_DT0d1_MESHwater0d15.obj_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_Trial06
./fast ./input_files/Ruehl2016_DT0d1_MESHwater0d15.obj_ELEMlinear_ALElinear_ALEPERIOD1_Trial06
./fast ./input_files/Ruehl2016_DT0d1_MESHwater0d15.obj_ELEMpseudo_quad_ALElinear_ALEPERIOD1_Trial06

./fast ./input_files/Ruehl2016_DT0d1_MESHwater0d15.obj_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1_Trial05
./fast ./input_files/Ruehl2016_DT0d1_MESHwater0d15.obj_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_Trial05
./fast ./input_files/Ruehl2016_DT0d1_MESHwater0d15.obj_ELEMlinear_ALElinear_ALEPERIOD1_Trial05
./fast ./input_files/Ruehl2016_DT0d1_MESHwater0d15.obj_ELEMpseudo_quad_ALElinear_ALEPERIOD1_Trial05

./fast ./input_files/Ruehl2016_DT0d1_MESHwater0d15.obj_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1_Trial04
./fast ./input_files/Ruehl2016_DT0d1_MESHwater0d15.obj_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_Trial04
./fast ./input_files/Ruehl2016_DT0d1_MESHwater0d15.obj_ELEMlinear_ALElinear_ALEPERIOD1_Trial04
./fast ./input_files/Ruehl2016_DT0d1_MESHwater0d15.obj_ELEMpseudo_quad_ALElinear_ALEPERIOD1_Trial04
