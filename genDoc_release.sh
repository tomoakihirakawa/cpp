#!/bin/bash

# python3 ./extract_comments_release.py README.md ./
# echo 'python3 ./extract_comments_release.py README.md ./'

(cd builds/build_sph; python3 ../../extract_comments_release.py README.md ./ ../../)
echo '(cd builds/build_sph; python3 ../../extract_comments_release.py README.md ./ ../../)'

(cd builds/build_bem; python3 ../../extract_comments_release.py README.md ./ ../../)
echo '(cd builds/build_bem; python3 ../../extract_comments_release.py README.md ./ ../../)'

(cd builds/build_ODE; python3 ../../extract_comments_release.py README.md ./ ../../)
echo '(cd builds/build_ODE; python3 ../../extract_comments_release.py README.md ./ ../../)'

(cd builds/build_divide_merge; python3 ../../extract_comments_release.py README.md ./ ../../)
echo '(cd builds/build_divide_merge; python3 ../../extract_comments_release.py README.md ./ ../../)'

(cd builds/build_spherical_harmonic; python3 ../../extract_comments_release.py README.md ./ ../../)
echo '(cd builds/build_spherical_harmonic; python3 ../../extract_comments_release.py README.md ./ ../../)'

(cd builds/build_system_of_linear_eqs; python3 ../../extract_comments_release.py README.md ./ ../../)
echo '(cd builds/build_system_of_linear_eqs; python3 ../../extract_comments_release.py README.md ./ ../../)'

(cd builds/build_root_finding; python3 ../../extract_comments_release.py README.md ./ ../../)
echo '(cd builds/build_root_finding; python3 ../../extract_comments_release.py README.md ./ ../../)'

(cd builds/build_pybind11; python3 ../../extract_comments_release.py README.md ./ ../../)
echo '(cd builds/build_pybind11; python3 ../../extract_comments_release.py README.md ./ ../../)'


python3 generate_readme.py