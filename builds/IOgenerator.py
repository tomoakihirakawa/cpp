'''DOC_EXTRACT 1_0_HOW_TO_MAKE_INPUT_FILES

# 入力ファイル生成 `input_generator.py`

This Python script generates input files for the BEM simulation codxxe. It supports various simulation cases and handles input file generation for each case.

## Usage

1. Make sure the required dependencies are installed.
2. Run the script using the following command:

```shell
$ python3 input_generator.py
```

Upon running the script, it will generate_input_files input files in JSON format for the specified simulation case. The input files are saved in the `./input_files/` directory.

## Customization

To customize the input file generation for a specific case, follow these steps:

1. Locate the `SimulationCase` variable in the script and set it to the desired case name, e.g., `"Kramer2021"`.
2. Add a new `case` block in the `match SimulationCase:` section to handle the new simulation case.
3. Define the required parameters for the simulation case within the new `case` block, following the examples provided in the script.
4. Update the `inputfiles` variable with the new input objects created for the custom case.

After customizing the script, run it again to generate_input_files the input files for the new case.

## Output

The script will generate_input_files input files in JSON format for the specified simulation case. The input files will be saved in the `./input_files/` directory. The generated input files can be used to run the BEM simulation.
'''


'''
プログラムを回す際に面倒な事は，入力ファイルの設定方法．
入力ファイルの作り方をドキュメントで示されても，具体的な例がないとわかりにくい．
例があっても，例と違う場合どうすればいいかなど，わからないことは多い．
このように，入力ファイルを生成するプログラムを作っておけば，その面倒をだいぶ解消できる．
'''

import copy
import platform
from math import pi
import json
import math
import os
from os.path import expanduser

home = expanduser("~")

white = '\033[90m'
red = '\033[91m'
yellow = '\033[93m'
blue = '\033[96m'
green = '\033[92m'
orange = '\033[33m'
magenta = '\033[95m'
coloroff = '\033[0m'
#boldcases
White = '\033[1;90m'
Red = '\033[1;91m'
Yellow = '\033[1;93m'
Blue = '\033[1;96m'
Green = '\033[1;92m'
Magenta = '\033[1;95m'


# output_directory and input_directory are automatically added to setting
def generate_input_files(inputfiles, setting, generate_in_out_directory, id):

    input_directory, output_directory = generate_in_out_directory(id)
    setting["output_directory"] = output_directory
    setting["input_files"] = [x["name"]+".json" for x in inputfiles]

    # if input_files is not found in setting, then add it
    if "input_files" not in setting:
        setting["input_files"] = [x["name"]+".json" for x in inputfiles]

    # @ -------------------------------------------------------- #
    # @           その他，water.json,tank.json などを出力           #
    # @ -------------------------------------------------------- #

    def does_file_exist(key, filename):
        if key == "objfile":
            if os.path.exists(filename) == False:
                return "❌" + Red + " file does not exist"+coloroff
            else:
                return "✅" + Green + " file exists"+coloroff
        return ''

    for INPUTS in inputfiles:
        NAMEJSON = INPUTS["name"]+'.json'
        print(blue,'-'*(40-len(NAMEJSON)), NAMEJSON, coloroff)
        for key, value in INPUTS.items():
            if key == "velocity":
                if isinstance(value, list) and "wave" in value[0]:  # Check if the value is a list and if "wave" is in the first element
                    print(f'{key: >{10}}', ':\t🌊', green, value, coloroff)
                elif value == "floating":
                    print(f'{key: >{10}}', ':\t🚢', green, value, coloroff)
                elif "flap" in value or "piston" in value:
                    print(f'{key: >{10}}', ':\t🌊', green, value, coloroff)
                else:
                    print(f'{key: >{10}}', ':\t', green, value, coloroff)
            elif value == "RigidBody":
                print(f'{key: >{10}}', ':\t🚀', magenta, value, coloroff)
            elif value == "Fluid":
                print(f'{key: >{10}}', ':\t💧', blue, value, coloroff)
            elif key == "type" and isinstance(value, (str, list)) and ("probe" in value or "sensor" in value or "gauge" in value):
                print(f'{key: >{10}}', ':\t📐', yellow, value, coloroff)
            elif key == "type" and isinstance(value, (str, list)) and ("absorb" in value or "damp"):
                print(f'{key: >{10}}', ':\t🏖️', orange, value, coloroff)
            elif key == "isFixed":
                if value == True:
                    print(f'{key: >{10}}', ':\t📌', orange, value, coloroff)
                else:
                    print(f'{key: >{10}}', ':\t', green, value, coloroff)
            elif "mooring" in key:
                print(f'{key: >{10}}', ':\t🔗', green, value, coloroff)
            else:
                print(f'{key: >{10}}', ':\t' + does_file_exist(key, value), white, value, coloroff)

        f = open(input_directory+"/"+INPUTS["name"]+".json", 'w')
        json.dump(INPUTS, f, ensure_ascii=True, indent=4)
        f.close()

    # @ -------------------------------------------------------- #
    # @                  setting.json を出力                      #
    # @ -------------------------------------------------------- #    

    print(Blue,'='*28,'setting.json', coloroff)
    for key, value in setting.items():
        print(f'{key: >{18}}', ':\t', green, value, coloroff)
    print(Blue,'-'*40, coloroff)
    f = open(input_directory+"/setting.json", 'w')
    json.dump(setting, f, ensure_ascii=True, indent=4)
    f.close()

    print("The directory for input files :", magenta, input_directory, coloroff)

# ---------------------------------------------------------------------------- #

try:
    import numpy as np
    numpy_available = True
except ImportError:
    numpy_available = False
    np = None

import numpy as np

def parse_vertex(line):
    parts = line.split()
    return np.array([float(parts[1]), float(parts[2]), float(parts[3])])

def read_obj(filename):
    vertices = []
    triangles = []
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('v '):
                vertex = parse_vertex(line)
                vertices.append(vertex)
            elif line.startswith('f '):
                parts = line.split()
                # Assuming the face is a triangle
                triangle = [int(parts[1].split('/')[0]), int(parts[2].split('/')[0]), int(parts[3].split('/')[0])]
                triangles.append(triangle)
    return len(vertices), len(triangles)