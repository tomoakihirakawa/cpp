#!/bin/sh

T_list=(0.75 0.74 0.73 0.72 0.71 0.7 0.69 0.68 0.67 0.66 0.65 0.64 0.63 0.62 0.61 0.6 0.59 0.58 0.57 0.56 0.55 0.54 0.53 0.52 0.51 0.5)

for t in ${T_list[@]}; do
    # python3 input_generator.py Tonegawa2024 -ALEPERIOD 1 -dt 0.025 -suffix B10d25 -outputdir /Volumes/home/BEM/Tonegawa2024_new
    python3.11 ./input_generator.py Tonegawa2024Experiment -T ${t} -dt 0.03 -s 0d25 -o /Volumes/home/BEM/Tonegawa2024Experiment
done

for t in ${T_list[@]}; do
    # python3 input_generator.py Tonegawa2024 -ALEPERIOD 1 -dt 0.025  -suffix B10d125 -outputdir /Volumes/home/BEM/Tonegawa2024_new
    python3.11 ./input_generator.py Tonegawa2024Experiment -T ${t} -dt 0.03 -s 0d125 -o /Volumes/home/BEM/Tonegawa2024Experiment
done

for t in ${T_list[@]}; do
    # python3 input_generator.py Tonegawa2024 -ALEPERIOD 1 -dt 0.025 -suffix B10d075 -outputdir /Volumes/home/BEM/Tonegawa2024_new
    python3.11 ./input_generator.py Tonegawa2024Experiment -T ${t} -dt 0.03 -s 0d075 -o /Volumes/home/BEM/Tonegawa2024Experiment
done

# ./main ./input_files/Tonegawa2024_T0d5_DT0d025_ELEMlinear_ALEquasi_quad_ALEPERIOD1_B10d075
# ./main ./input_files/Tonegawa2024_T0d51_DT0d025_ELEMlinear_ALEquasi_quad_ALEPERIOD1_B10d075
# ./main ./input_files/Tonegawa2024_T0d52_DT0d025_ELEMlinear_ALEquasi_quad_ALEPERIOD1_B10d075
# ./main ./input_files/Tonegawa2024_T0d53_DT0d025_ELEMlinear_ALEquasi_quad_ALEPERIOD1_B10d075
# ./main ./input_files/Tonegawa2024_T0d54_DT0d025_ELEMlinear_ALEquasi_quad_ALEPERIOD1_B10d075

# ./main ./input_files/Tonegawa2024_T0d55_DT0d025_ELEMlinear_ALEquasi_quad_ALEPERIOD1_B10d075
# ./main ./input_files/Tonegawa2024_T0d56_DT0d025_ELEMlinear_ALEquasi_quad_ALEPERIOD1_B10d075
# ./main ./input_files/Tonegawa2024_T0d57_DT0d025_ELEMlinear_ALEquasi_quad_ALEPERIOD1_B10d075
# ./main ./input_files/Tonegawa2024_T0d58_DT0d025_ELEMlinear_ALEquasi_quad_ALEPERIOD1_B10d075
# ./main ./input_files/Tonegawa2024_T0d59_DT0d025_ELEMlinear_ALEquasi_quad_ALEPERIOD1_B10d075

# ./main ./input_files/Tonegawa2024_T0d6_DT0d025_ELEMlinear_ALEquasi_quad_ALEPERIOD1_B10d075
# ./main ./input_files/Tonegawa2024_T0d61_DT0d025_ELEMlinear_ALEquasi_quad_ALEPERIOD1_B10d075
# ./main ./input_files/Tonegawa2024_T0d62_DT0d025_ELEMlinear_ALEquasi_quad_ALEPERIOD1_B10d075
# ./main ./input_files/Tonegawa2024_T0d63_DT0d025_ELEMlinear_ALEquasi_quad_ALEPERIOD1_B10d075
# ./main ./input_files/Tonegawa2024_T0d64_DT0d025_ELEMlinear_ALEquasi_quad_ALEPERIOD1_B10d075
# ./main ./input_files/Tonegawa2024_T0d65_DT0d025_ELEMlinear_ALEquasi_quad_ALEPERIOD1_B10d075

# ./main ./input_files/Tonegawa2024_T0d66_DT0d025_ELEMlinear_ALEquasi_quad_ALEPERIOD1_B10d075
# ./main ./input_files/Tonegawa2024_T0d67_DT0d025_ELEMlinear_ALEquasi_quad_ALEPERIOD1_B10d075
# ./main ./input_files/Tonegawa2024_T0d68_DT0d025_ELEMlinear_ALEquasi_quad_ALEPERIOD1_B10d075
# ./main ./input_files/Tonegawa2024_T0d69_DT0d025_ELEMlinear_ALEquasi_quad_ALEPERIOD1_B10d075
# ./main ./input_files/Tonegawa2024_T0d7_DT0d025_ELEMlinear_ALEquasi_quad_ALEPERIOD1_B10d075

# ./fast ./input_files/Tonegawa2024Experiment_T0d50_DT0d03_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d075
# ./fast ./input_files/Tonegawa2024Experiment_T0d52_DT0d03_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d075
# ./fast ./input_files/Tonegawa2024Experiment_T0d54_DT0d03_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d075
# ./fast ./input_files/Tonegawa2024Experiment_T0d56_DT0d03_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d075
# ./fast ./input_files/Tonegawa2024Experiment_T0d58_DT0d03_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d075
# ./fast ./input_files/Tonegawa2024Experiment_T0d60_DT0d03_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d075
# ./fast ./input_files/Tonegawa2024Experiment_T0d62_DT0d03_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d075
# ./fast ./input_files/Tonegawa2024Experiment_T0d64_DT0d03_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d075
# ./fast ./input_files/Tonegawa2024Experiment_T0d66_DT0d03_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d075
# ./fast ./input_files/Tonegawa2024Experiment_T0d68_DT0d03_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d075
# ./fast ./input_files/Tonegawa2024Experiment_T0d70_DT0d03_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d075
# ./fast ./input_files/Tonegawa2024Experiment_T0d72_DT0d03_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d075
# ./fast ./input_files/Tonegawa2024Experiment_T0d74_DT0d03_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d075
# ./fast ./input_files/Tonegawa2024Experiment_T0d76_DT0d03_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d075

# ./fast ./input_files/Tonegawa2024Experiment_T0d50_DT0d03_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d125
# ./fast ./input_files/Tonegawa2024Experiment_T0d52_DT0d03_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d125
# ./fast ./input_files/Tonegawa2024Experiment_T0d54_DT0d03_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d125
# ./fast ./input_files/Tonegawa2024Experiment_T0d56_DT0d03_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d125
# ./fast ./input_files/Tonegawa2024Experiment_T0d58_DT0d03_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d125
# ./fast ./input_files/Tonegawa2024Experiment_T0d60_DT0d03_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d125
# ./fast ./input_files/Tonegawa2024Experiment_T0d62_DT0d03_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d125
# ./fast ./input_files/Tonegawa2024Experiment_T0d64_DT0d03_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d125
# ./fast ./input_files/Tonegawa2024Experiment_T0d66_DT0d03_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d125
# ./fast ./input_files/Tonegawa2024Experiment_T0d68_DT0d03_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d125
# ./fast ./input_files/Tonegawa2024Experiment_T0d70_DT0d03_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d125
# ./fast ./input_files/Tonegawa2024Experiment_T0d72_DT0d03_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d125
# ./fast ./input_files/Tonegawa2024Experiment_T0d74_DT0d03_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d125
# ./fast ./input_files/Tonegawa2024Experiment_T0d76_DT0d03_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_0d125