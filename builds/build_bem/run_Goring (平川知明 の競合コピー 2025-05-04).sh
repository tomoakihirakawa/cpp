#!/bin/bash

# python3.11 input_generator.py Goring1979 -m water0d06refined.obj -dt 0.05 -e pseudo_quad -ALE pseudo_quad -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501
# python3.11 input_generator.py Goring1979 -m water0d05refined.obj -dt 0.05 -e pseudo_quad -ALE pseudo_quad -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501

# python3.11 input_generator.py Goring1979 -m water0d06refined.obj -dt 0.05 -e pseudo_quad -ALE linear -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501
# python3.11 input_generator.py Goring1979 -m water0d05refined.obj -dt 0.05 -e pseudo_quad -ALE linear -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501

# python3.11 input_generator.py Goring1979 -m water0d06refined.obj -dt 0.05 -e linear -ALE pseudo_quad -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501
# python3.11 input_generator.py Goring1979 -m water0d05refined.obj -dt 0.05 -e linear -ALE pseudo_quad -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501

# python3.11 input_generator.py Goring1979 -m water0d06refined.obj -dt 0.05 -e linear -ALE linear -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501
# python3.11 input_generator.py Goring1979 -m water0d05refined.obj -dt 0.05 -e linear -ALE linear -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501

# ./fast ./input_files/Goring1979_DT0d05_MESHwater0d06refined.obj_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1
# ./fast ./input_files/Goring1979_DT0d05_MESHwater0d05refined.obj_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1
# ./fast ./input_files/Goring1979_DT0d05_MESHwater0d06refined.obj_ELEMpseudo_quad_ALElinear_ALEPERIOD1
# ./fast ./input_files/Goring1979_DT0d05_MESHwater0d05refined.obj_ELEMpseudo_quad_ALElinear_ALEPERIOD1

# ./fast ./input_files/Goring1979_DT0d05_MESHwater0d06refined.obj_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
# ./fast ./input_files/Goring1979_DT0d05_MESHwater0d05refined.obj_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
# ./fast ./input_files/Goring1979_DT0d05_MESHwater0d06refined.obj_ELEMlinear_ALElinear_ALEPERIOD1
# ./fast ./input_files/Goring1979_DT0d05_MESHwater0d05refined.obj_ELEMlinear_ALElinear_ALEPERIOD1

# python3.11 input_generator.py Ruehl2016 -m water0d2.obj -dt 0.1 -e pseudo_quad -ALE pseudo_quad -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501
# python3.11 input_generator.py Ruehl2016 -m water0d2.obj -dt 0.1 -e pseudo_quad -ALE linear -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501
# python3.11 input_generator.py Ruehl2016 -m water0d2.obj -dt 0.1 -e linear -ALE pseudo_quad -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501
# python3.11 input_generator.py Ruehl2016 -m water0d2.obj -dt 0.1 -e linear -ALE linear -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501

# ./fast ./input_files/Ruehl2016_DT0d1_MESHwater0d2.obj_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
# ./fast ./input_files/Ruehl2016_DT0d1_MESHwater0d2.obj_ELEMlinear_ALElinear_ALEPERIOD1
# ./fast ./input_files/Ruehl2016_DT0d1_MESHwater0d2.obj_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1
# ./fast ./input_files/Ruehl2016_DT0d1_MESHwater0d2.obj_ELEMpseudo_quad_ALElinear_ALEPERIOD1

python3.11 input_generator.py Ruehl2016 -m water0d075.obj -dt 0.05 -e pseudo_quad -ALE pseudo_quad -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501
python3.11 input_generator.py Ruehl2016 -m water0d075.obj -dt 0.05 -e pseudo_quad -ALE linear -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501
python3.11 input_generator.py Ruehl2016 -m water0d075.obj -dt 0.05 -e linear -ALE pseudo_quad -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501
python3.11 input_generator.py Ruehl2016 -m water0d075.obj -dt 0.05 -e linear -ALE linear -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501

python3.11 input_generator.py Ruehl2016 -m water0d1.obj -dt 0.05 -e pseudo_quad -ALE pseudo_quad -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501
python3.11 input_generator.py Ruehl2016 -m water0d1.obj -dt 0.05 -e pseudo_quad -ALE linear -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501
python3.11 input_generator.py Ruehl2016 -m water0d1.obj -dt 0.05 -e linear -ALE pseudo_quad -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501
python3.11 input_generator.py Ruehl2016 -m water0d1.obj -dt 0.05 -e linear -ALE linear -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501

./fast ./input_files/Ruehl2016_DT0d05_MESHwater0d075.obj_ELEMlinear_ALElinear_ALEPERIOD1
./fast ./input_files/Ruehl2016_DT0d05_MESHwater0d075.obj_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
./fast ./input_files/Ruehl2016_DT0d05_MESHwater0d075.obj_ELEMpseudo_quad_ALElinear_ALEPERIOD1
./fast ./input_files/Ruehl2016_DT0d05_MESHwater0d075.obj_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1

./fast ./input_files/Ruehl2016_DT0d05_MESHwater0d1.obj_ELEMlinear_ALElinear_ALEPERIOD1
./fast ./input_files/Ruehl2016_DT0d05_MESHwater0d1.obj_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
./fast ./input_files/Ruehl2016_DT0d05_MESHwater0d1.obj_ELEMpseudo_quad_ALElinear_ALEPERIOD1
./fast ./input_files/Ruehl2016_DT0d05_MESHwater0d1.obj_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1

# python3.11 input_generator.py Ruehl2016 -m water0d125.obj -dt 0.05 -e pseudo_quad -ALE pseudo_quad -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501
# python3.11 input_generator.py Ruehl2016 -m water0d125.obj -dt 0.05 -e pseudo_quad -ALE linear -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501
# python3.11 input_generator.py Ruehl2016 -m water0d125.obj -dt 0.05 -e linear -ALE pseudo_quad -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501
# python3.11 input_generator.py Ruehl2016 -m water0d125.obj -dt 0.05 -e linear -ALE linear -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250501

# ./fast ./input_files/Ruehl2016_DT0d05_MESHwater0d125.obj_ELEMlinear_ALElinear_ALEPERIOD1
# ./fast ./input_files/Ruehl2016_DT0d05_MESHwater0d125.obj_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
# ./fast ./input_files/Ruehl2016_DT0d05_MESHwater0d125.obj_ELEMpseudo_quad_ALElinear_ALEPERIOD1
# ./fast ./input_files/Ruehl2016_DT0d05_MESHwater0d125.obj_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1

# ./fast ./input_files/Ruehl2016_DT0d05_MESHwater0d1.obj_ELEMlinear_ALElinear_ALEPERIOD1
# ./fast ./input_files/Ruehl2016_DT0d05_MESHwater0d1.obj_ELEMlinear_ALEpseudo_quad_ALEPERIOD1

# python3.11 input_generator.py Goring1979 -m water0d09refined.obj -dt 0.03 -e pseudo_quad -ALE pseudo_quad -s 20250329 -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250329
# python3.11 input_generator.py Goring1979 -m water0d09refined.obj -dt 0.03 -e linear -ALE pseudo_quad -s 20250329 -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250329
# python3.11 input_generator.py Goring1979 -m water0d09refined.obj -dt 0.03 -e pseudo_quad -ALE linear -s 20250329 -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250329
# python3.11 input_generator.py Goring1979 -m water0d09refined.obj -dt 0.03 -e linear -ALE linear -s 20250329 -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250329

# for dt in 0.03 0.06; do
#    python3.11 input_generator.py Goring1979 -m water0d06refined.obj -dt $dt -e pseudo_quad -ALE pseudo_quad -s 20250329 -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250329
#    python3.11 input_generator.py Goring1979 -m water0d06refined.obj -dt $dt -e linear -ALE pseudo_quad -s 20250329 -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250329
#    python3.11 input_generator.py Goring1979 -m water0d06refined.obj -dt $dt -e pseudo_quad -ALE linear -s 20250329 -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250329
#    python3.11 input_generator.py Goring1979 -m water0d06refined.obj -dt $dt -e linear -ALE linear -s 20250329 -ALEPERIOD 1 -o /Volumes/home/BEM/benchmark20250329
# done

# # やはり造波性能の問題がるのかもしれない．Grilli2004も造波された波の形は実験とよく一致していない．

# ./fast ./input_files/Goring1979_DT0d03_MESHwater0d09refined.obj_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1_20250329
# ./fast ./input_files/Goring1979_DT0d03_MESHwater0d09refined.obj_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_20250329
# ./fast ./input_files/Goring1979_DT0d03_MESHwater0d09refined.obj_ELEMlinear_ALElinear_ALEPERIOD1_20250329
# ./fast ./input_files/Goring1979_DT0d03_MESHwater0d09refined.obj_ELEMpseudo_quad_ALElinear_ALEPERIOD1_20250329

# ./fast ./input_files/Goring1979_DT0d03_MESHwater0d06refined.obj_ELEMpseudo_quad_ALEpseudo_quad_ALEPERIOD1_20250329
# ./fast ./input_files/Goring1979_DT0d03_MESHwater0d06refined.obj_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_20250329
# ./fast ./input_files/Goring1979_DT0d03_MESHwater0d06refined.obj_ELEMlinear_ALElinear_ALEPERIOD1_20250329
# ./fast ./input_files/Goring1979_DT0d03_MESHwater0d06refined.obj_ELEMpseudo_quad_ALElinear_ALEPERIOD1_20250329
