#!/bin/sh

T_list=$(seq 5 0.5 9)

# python3.11 input_generator.py Tonegawa2024Akita -dt 0.2 -T 5 -H 2 -o /Volumes/home/BEM/Tonegawa2024Akita/ -s 3MW_MP30

for t in ${T_list[@]}; do
   python3.11 input_generator.py Tonegawa2024Akita -dt 0.4 -T ${t} -H 1 -s 3MW_MP30_new -o /Volumes/home/BEM/Tonegawa2024Akita/
   python3.11 input_generator.py Tonegawa2024Akita -dt 0.4 -T ${t} -H 1 -s 3MW_MP15_new -o /Volumes/home/BEM/Tonegawa2024Akita/
   python3.11 input_generator.py Tonegawa2024Akita -dt 0.4 -T ${t} -H 1 -s 3MW_MP9_new -o /Volumes/home/BEM/Tonegawa2024Akita/

   python3.11 input_generator.py Tonegawa2024Akita -dt 0.4 -T ${t} -H 1 -s 10MW_MP30_new -o /Volumes/home/BEM/Tonegawa2024Akita/
   python3.11 input_generator.py Tonegawa2024Akita -dt 0.4 -T ${t} -H 1 -s 10MW_MP15_new -o /Volumes/home/BEM/Tonegawa2024Akita/
   python3.11 input_generator.py Tonegawa2024Akita -dt 0.4 -T ${t} -H 1 -s 10MW_MP9_new -o /Volumes/home/BEM/Tonegawa2024Akita/
done

# ./fast ./input_files/Tonegawa2024Akita_H1d0_T5d0_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_10MW_MP30_new
# ./fast ./input_files/Tonegawa2024Akita_H1d0_T5d5_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_10MW_MP30_new
# ./fast ./input_files/Tonegawa2024Akita_H1d0_T6d0_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_10MW_MP30_new
# ./fast ./input_files/Tonegawa2024Akita_H1d0_T6d5_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_10MW_MP30_new
# ./fast ./input_files/Tonegawa2024Akita_H1d0_T7d0_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_10MW_MP30_new
# ./fast ./input_files/Tonegawa2024Akita_H1d0_T7d5_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_10MW_MP30_new
# ./fast ./input_files/Tonegawa2024Akita_H1d0_T8d0_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_10MW_MP30_new
# ./fast ./input_files/Tonegawa2024Akita_H1d0_T8d5_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_10MW_MP30_new
# ./fast ./input_files/Tonegawa2024Akita_H1d0_T9d0_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_10MW_MP30_new

# ./fast ./input_files/Tonegawa2024Akita_H1d0_T5d0_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_3MW_MP15_new
# ./fast ./input_files/Tonegawa2024Akita_H1d0_T5d5_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_3MW_MP15_new
# ./fast ./input_files/Tonegawa2024Akita_H1d0_T6d0_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_3MW_MP15_new
# ./fast ./input_files/Tonegawa2024Akita_H1d0_T6d5_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_3MW_MP15_new
# ./fast ./input_files/Tonegawa2024Akita_H1d0_T7d0_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_3MW_MP15_new
# ./fast ./input_files/Tonegawa2024Akita_H1d0_T7d5_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_3MW_MP15_new
# ./fast ./input_files/Tonegawa2024Akita_H1d0_T8d0_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_3MW_MP15_new
# ./fast ./input_files/Tonegawa2024Akita_H1d0_T8d5_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_3MW_MP15_new
# ./fast ./input_files/Tonegawa2024Akita_H1d0_T9d0_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_3MW_MP15_new

# ./fast ./input_files/Tonegawa2024Akita_H1d0_T5d0_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_3MW_MP9
# ./fast ./input_files/Tonegawa2024Akita_H1d0_T5d2_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_3MW_MP9
# ./fast ./input_files/Tonegawa2024Akita_H1d0_T5d4_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_3MW_MP9
# ./fast ./input_files/Tonegawa2024Akita_H1d0_T5d6_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_3MW_MP9

# ./fast ./input_files/Tonegawa2024Akita_H1d0_T5d8_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_3MW_MP9
# ./fast ./input_files/Tonegawa2024Akita_H1d0_T6d0_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_3MW_MP9
# ./fast ./input_files/Tonegawa2024Akita_H1d0_T6d2_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_3MW_MP9
# ./fast ./input_files/Tonegawa2024Akita_H1d0_T6d4_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_3MW_MP9

# ./fast ./input_files/Tonegawa2024Akita_H1d0_T6d6_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_3MW_MP9
# ./fast ./input_files/Tonegawa2024Akita_H1d0_T6d8_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_3MW_MP9
# ./fast ./input_files/Tonegawa2024Akita_H1d0_T7d0_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_3MW_MP9
# ./fast ./input_files/Tonegawa2024Akita_H1d0_T7d2_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_3MW_MP9

# ./fast ./input_files/Tonegawa2024Akita_H1d0_T7d4_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_3MW_MP9
# ./fast ./input_files/Tonegawa2024Akita_H1d0_T7d6_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_3MW_MP9
# ./fast ./input_files/Tonegawa2024Akita_H1d0_T7d8_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_3MW_MP9
# ./fast ./input_files/Tonegawa2024Akita_H1d0_T8d0_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_3MW_MP9

# ./fast ./input_files/Tonegawa2024Akita_H1d0_T8d2_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_10MW_MP9
# ./fast ./input_files/Tonegawa2024Akita_H1d0_T8d4_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_10MW_MP9
# ./fast ./input_files/Tonegawa2024Akita_H1d0_T8d6_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_10MW_MP9
# ./fast ./input_files/Tonegawa2024Akita_H1d0_T8d8_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_10MW_MP9
# ./fast ./input_files/Tonegawa2024Akita_H1d0_T9d0_DT0d4_ELEMlinear_ALEpseudo_quad_ALEPERIOD1_10MW_MP9