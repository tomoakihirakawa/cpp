#!/bin/sh

# outputdir=/Volumes/home/BEM/benchmark202404

# set home directory
outputdir=/Users/tomoaki

# case = ${case}
case=Tanizawa1996
# case=Palm2016

for H in 0.05 0.3; do
    # H as number
    # 擬2次要素と線形要素の比較 ALEは線形要素上で行う　メッシュも細かくして収束を確認

    python3.11 input_generator.py -case ${case} -mesh water_no_float0d08 -element pseudo_quad -wavemaker flap -dt 0.03 -H ${H} -ALE linear -outputdir ${outputdir}
    python3.11 input_generator.py -case ${case} -mesh water_no_float0d08 -element linear -wavemaker flap -dt 0.03 -H ${H} -ALE linear -outputdir ${outputdir}
    python3.11 input_generator.py -case ${case} -mesh water_no_float0d07 -element pseudo_quad -wavemaker flap -dt 0.03 -H ${H} -ALE linear -outputdir ${outputdir}
    python3.11 input_generator.py -case ${case} -mesh water_no_float0d07 -element linear -wavemaker flap -dt 0.03 -H ${H} -ALE linear -outputdir ${outputdir}

    # 擬2次要素と線形要素の比較 ALEは擬2次要素上で行う　メッシュも細かくして収束を確認
    python3.11 input_generator.py -case ${case} -mesh water_no_float0d08 -element pseudo_quad -wavemaker flap -dt 0.03 -H ${H} -ALE pseudo_quad -outputdir ${outputdir}
    python3.11 input_generator.py -case ${case} -mesh water_no_float0d08 -element linear -wavemaker flap -dt 0.03 -H ${H} -ALE pseudo_quad -outputdir ${outputdir}
    python3.11 input_generator.py -case ${case} -mesh water_no_float0d07 -element pseudo_quad -wavemaker flap -dt 0.03 -H ${H} -ALE pseudo_quad -outputdir ${outputdir}
    python3.11 input_generator.py -case ${case} -mesh water_no_float0d07 -element linear -wavemaker flap -dt 0.03 -H ${H} -ALE pseudo_quad -outputdir ${outputdir}

    # タイムステップを短くして収束を確認
    python3.11 input_generator.py -case ${case} -mesh water_no_float0d08 -element pseudo_quad -wavemaker flap -dt 0.05 -H ${H} -ALE linear -outputdir ${outputdir}
    python3.11 input_generator.py -case ${case} -mesh water_no_float0d08 -element linear -wavemaker flap -dt 0.05 -H ${H} -ALE linear -outputdir ${outputdir}
    python3.11 input_generator.py -case ${case} -mesh water_no_float0d07 -element pseudo_quad -wavemaker flap -dt 0.05 -H ${H} -ALE linear -outputdir ${outputdir}
    python3.11 input_generator.py -case ${case} -mesh water_no_float0d07 -element linear -wavemaker flap -dt 0.05 -H ${H} -ALE linear -outputdir ${outputdir}

    # タイムステップを短くして収束を確認
    python3.11 input_generator.py -case ${case} -mesh water_no_float0d08 -element pseudo_quad -wavemaker flap -dt 0.05 -H ${H} -ALE pseudo_quad -outputdir ${outputdir}
    python3.11 input_generator.py -case ${case} -mesh water_no_float0d08 -element linear -wavemaker flap -dt 0.05 -H ${H} -ALE pseudo_quad -outputdir ${outputdir}
    python3.11 input_generator.py -case ${case} -mesh water_no_float0d07 -element pseudo_quad -wavemaker flap -dt 0.05 -H ${H} -ALE pseudo_quad -outputdir ${outputdir}
    python3.11 input_generator.py -case ${case} -mesh water_no_float0d07 -element linear -wavemaker flap -dt 0.05 -H ${H} -ALE pseudo_quad -outputdir ${outputdir}

done

# ./main ./input_files/${case}_H0d3_L1d8_WAVEflap_MESHwater_no_float0d08_DT0d05_ELEMpseudo_quad_ALEpseudo_quad

# ----------------------------------- H005 ----------------------------------- #

# mac studio
# ./main ./input_files/${case}_H0d05_L1d8_WAVEflap_MESHwater_no_float0d08_DT0d05_ELEMpseudo_quad_ALEpseudo_quad
# ./main ./input_files/${case}_H0d05_L1d8_WAVEflap_MESHwater_no_float0d08_DT0d05_ELEMlinear_ALEpseudo_quad
# ./main ./input_files/${case}_H0d05_L1d8_WAVEflap_MESHwater_no_float0d08_DT0d05_ELEMpseudo_quad_ALElinear
# ./main ./input_files/${case}_H0d05_L1d8_WAVEflap_MESHwater_no_float0d08_DT0d05_ELEMlinear_ALElinear

# ./main ./input_files/${case}_H0d05_L1d8_WAVEflap_MESHwater_no_float0d08_DT0d03_ELEMpseudo_quad_ALEpseudo_quad # prandtl 142
# ./main ./input_files/${case}_H0d05_L1d8_WAVEflap_MESHwater_no_float0d08_DT0d03_ELEMlinear_ALEpseudo_quad # Daiki 140
# ./main ./input_files/${case}_H0d05_L1d8_WAVEflap_MESHwater_no_float0d08_DT0d03_ELEMpseudo_quad_ALElinear #macbook 
# ./main ./input_files/${case}_H0d05_L1d8_WAVEflap_MESHwater_no_float0d08_DT0d03_ELEMlinear_ALElinear #ishiwaka 141

# ------------------------------------ H01 ----------------------------------- #

# ./main ./input_files/${case}_H0d1_L1d8_WAVEflap_MESHwater_no_float0d08_DT0d05_ELEMpseudo_quad_ALEpseudo_quad
# ./main ./input_files/${case}_H0d1_L1d8_WAVEflap_MESHwater_no_float0d08_DT0d05_ELEMlinear_ALEpseudo_quad
# ./main ./input_files/${case}_H0d1_L1d8_WAVEflap_MESHwater_no_float0d08_DT0d05_ELEMpseudo_quad_ALElinear
# ./main ./input_files/${case}_H0d1_L1d8_WAVEflap_MESHwater_no_float0d08_DT0d05_ELEMlinear_ALElinear

# ./main ./input_files/${case}_H0d1_L1d8_WAVEflap_MESHwater_no_float0d08_DT0d03_ELEMpseudo_quad_ALEpseudo_quad # prandtl 142
# ./main ./input_files/${case}_H0d1_L1d8_WAVEflap_MESHwater_no_float0d08_DT0d03_ELEMlinear_ALEpseudo_quad# Daiki 140
# ./main ./input_files/${case}_H0d1_L1d8_WAVEflap_MESHwater_no_float0d08_DT0d03_ELEMpseudo_quad_ALElinear #macbook 
# ./main ./input_files/${case}_H0d1_L1d8_WAVEflap_MESHwater_no_float0d08_DT0d03_ELEMlinear_ALElinear #ishiwaka 141
