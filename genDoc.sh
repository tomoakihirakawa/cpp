#!/bin/bash

python3 ./extract_comments.py README.md ./
echo 'python3 ./extract_comments.py README.md ./'

(cd builds/build_sph; python3 ../../extract_comments.py README.md ./ ../../)
echo '(cd builds/build_sph; python3 ../../extract_comments.py README.md ./ ../../)'

(cd builds/build_bem; python3 ../../extract_comments.py README.md ./ ../../)
echo '(cd builds/build_bem; python3 ../../extract_comments.py README.md ./ ../../)'

(cd builds/build_ODE; python3 ../../extract_comments.py README.md ./ ../../)
echo '(cd builds/build_ODE; python3 ../../extract_comments.py README.md ./ ../../)'

(cd builds/build_divide_merge; python3 ../../extract_comments.py README.md ./ ../../)
echo '(cd builds/build_divide_merge; python3 ../../extract_comments.py README.md ./ ../../)'