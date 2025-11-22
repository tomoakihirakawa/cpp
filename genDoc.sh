#!/bin/bash

# python3.11 ./extract_comments.py README.md ./
# echo 'python3.11 ./extract_comments.py README.md ./'

cd /Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/

python3.11 bib2json.py

# SPHのREADME.mdを生成
(
   cd builds/build_sph
   python3.11 ../../extract_comments.py README.md -main ./ -include ../../include/
)
# 学生用のREADME.mdを生成
# (
#    cd builds/build_sph
#    python3.11 ../../extract_comments.py README_FOR_STUDENTS.md -key EXTRACT_README_FOR_STUDENTS -main ./
# )

# 境界要素法のREADME.mdを生成
(
   cd builds/build_bem
   python3.11 ../../extract_comments.py README.md -main ./ -include ../../include/
)

# 学生用のREADME.mdを生成
# (
#    cd builds/build_bem
#    python3.11 ../../extract_comments.py README_FOR_STUDENTS.md -key EXTRACT_README_FOR_STUDENTS -main ./
# )
# (
#    cd builds/build_bem
#    python3.11 ../../extract_comments.py REVIEW_NOTE0.md -key EXTRACT_REVIEW_NOTE0 -main ./
# )
# (
#    cd builds/build_bem
#    python3.11 ../../extract_comments.py REVIEW_NOTE1.md -key EXTRACT_REVIEW_NOTE1 -main ./
# )

(
   cd builds/build_ODE
   python3.11 ../../extract_comments.py README.md -main ./ -include ../../include/
)

(
   cd builds/build_remesh
   python3.11 ../../extract_comments.py README.md -main ./ -include ../../include/
)

(
   cd builds/build_spherical_harmonic
   python3.11 ../../extract_comments.py README.md -main ./ -include ../../include/
)

(
   cd builds/build_system_of_linear_eqs
   python3.11 ../../extract_comments.py README.md -main ./ -include ../../include/
)

(
   cd builds/build_root_finding
   python3.11 ../../extract_comments.py README.md -main ./ -include ../../include/
)

(
   cd builds/build_pybind11
   python3.11 ../../extract_comments.py README.md -main ./ -include ../../include/
)

(
   cd builds/build_interpolation
   python3.11 ../../extract_comments.py README.md -main ./ -include ../../include/
)

(
   cd builds/build_quaternion
   python3.11 ../../extract_comments.py README.md -main ./ -include ../../include/
)

(
   cd builds/build_Network
   python3.11 ../../extract_comments.py README.md -main ./ -include ../../include/
)

(
   cd builds/build_JSON
   python3.11 ../../extract_comments.py README.md -main ./ -include ../../include/
)

(
   cd builds/build_eigen_value
   python3.11 ../../extract_comments.py README.md -main ./ -include ../../include/
)

(
   cd builds/build_integration
   python3.11 ../../extract_comments.py README.md -main ./ -include ../../include/
)

(
   cd builds/build_cable
   python3.11 ../../extract_comments.py README.md -main ./ -include ../../include/
)

(
   cd builds/build_Fourier
   python3.11 ../../extract_comments.py README.md -main ./ -include ../../include/
)

(
   cd builds/build_tensors
   python3.11 ../../extract_comments.py README.md -main ./ -include ../../include/
)

python3.11 generate_readme.py

cd -

# ---------------------------------------------------------------------------- #
