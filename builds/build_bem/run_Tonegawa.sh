#!/bin/sh

# T_list=(0.75 0.74 0.73 0.72 0.71 0.7 0.69 0.68 0.67 0.66 0.65 0.64 0.63 0.62 0.61 0.6 0.59 0.58 0.57 0.56 0.55 0.54 0.53 0.52 0.51 0.5)
T_list=(0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75)

for t in ${T_list[@]}; do
    # python3 input_generator.py Tonegawa2024 -ALEPERIOD 1 -dt 0.025 -suffix B10d25 -outputdir /Volumes/home/BEM/Tonegawa2024_new
    python3.11 ./input_generator.py Tonegawa2024Experiment -T ${t} -dt 0.05 -s 0d25 -o /Volumes/home/BEM/Tonegawa2024Experiment
done

for t in ${T_list[@]}; do
    # python3 input_generator.py Tonegawa2024 -ALEPERIOD 1 -dt 0.025  -suffix B10d125 -outputdir /Volumes/home/BEM/Tonegawa2024_new
    python3.11 ./input_generator.py Tonegawa2024Experiment -T ${t} -dt 0.05 -s 0d125 -o /Volumes/home/BEM/Tonegawa2024Experiment
done

for t in ${T_list[@]}; do
    # python3 input_generator.py Tonegawa2024 -ALEPERIOD 1 -dt 0.025 -suffix B10d075 -outputdir /Volumes/home/BEM/Tonegawa2024_new
    python3.11 ./input_generator.py Tonegawa2024Experiment -T ${t} -dt 0.05 -s 0d075 -o /Volumes/home/BEM/Tonegawa2024Experiment
done

# ./fast ./input_files/Tonegawa2024Experiment_T0d30_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d075
# ./fast ./input_files/Tonegawa2024Experiment_T0d35_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d075
# ./fast ./input_files/Tonegawa2024Experiment_T0d40_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d075
# ./fast ./input_files/Tonegawa2024Experiment_T0d45_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d075
# ./fast ./input_files/Tonegawa2024Experiment_T0d50_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d075

# ./fast ./input_files/Tonegawa2024Experiment_T0d55_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d075
# ./fast ./input_files/Tonegawa2024Experiment_T0d60_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d075
# ./fast ./input_files/Tonegawa2024Experiment_T0d65_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d075
# ./fast ./input_files/Tonegawa2024Experiment_T0d70_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d075
# ./fast ./input_files/Tonegawa2024Experiment_T0d75_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d075

# ./fast ./input_files/Tonegawa2024Experiment_T0d30_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d125
# ./fast ./input_files/Tonegawa2024Experiment_T0d35_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d125
# ./fast ./input_files/Tonegawa2024Experiment_T0d40_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d125
# ./fast ./input_files/Tonegawa2024Experiment_T0d45_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d125
# ./fast ./input_files/Tonegawa2024Experiment_T0d50_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d125

# ./fast ./input_files/Tonegawa2024Experiment_T0d55_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d125
# ./fast ./input_files/Tonegawa2024Experiment_T0d60_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d125
# ./fast ./input_files/Tonegawa2024Experiment_T0d65_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d125
# ./fast ./input_files/Tonegawa2024Experiment_T0d70_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d125
# ./fast ./input_files/Tonegawa2024Experiment_T0d75_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d125

# ./fast ./input_files/Tonegawa2024Experiment_T0d30_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d125
# ./fast ./input_files/Tonegawa2024Experiment_T0d35_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d125
# ./fast ./input_files/Tonegawa2024Experiment_T0d40_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d125
# ./fast ./input_files/Tonegawa2024Experiment_T0d45_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d125
# ./fast ./input_files/Tonegawa2024Experiment_T0d50_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d125

# ./fast ./input_files/Tonegawa2024Experiment_T0d55_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d25
# ./fast ./input_files/Tonegawa2024Experiment_T0d60_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d25
# ./fast ./input_files/Tonegawa2024Experiment_T0d65_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d25
# ./fast ./input_files/Tonegawa2024Experiment_T0d70_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d25
# ./fast ./input_files/Tonegawa2024Experiment_T0d75_DT0d05_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d25