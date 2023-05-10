from os.path import expanduser
import os
import math
import json
from math import pi
import platform

home = expanduser("~")

if platform.system() == "Linux":
    program_home = home + "/research"
else:
    program_home = home + "/Dropbox/code"

rho = 1000.
g = 9.81

'''
プログラムを回す際に面倒な事は，入力ファイルの設定方法．
入力ファイルの作り方をドキュメントで示されても，具体的な例がないとわかりにくい．
例があっても，例と違う場合どうすればいいかなど，わからないことは多い．
このように，入力ファイルを生成するプログラムを作っておけば，その面倒をだいぶ解消できる．
'''

# ---------------------------------------------------------------------------- #

SimulationCase = "Kamra2019"
id = ""
match SimulationCase:
    case "static_pressure":

        objfolder = program_home + "/cpp/obj/2023Tanabe/test20221219"

        water = {"name": "water",
                 "type": "Fluid",
                 "output_vtu_file_name": "water",  # 拡張子はいらない
                 "output_pvd_file_name": "water",  # 拡張子はいらない
                 "objfile": objfolder + "/water.obj"
                 }

        wavetank = {"name": "wavetank",
                    "type": "RigidBody",
                    "output_vtu_file_name": "wavetank",  # 拡張子はいらない
                    "output_pvd_file_name": "wavetank",  # 拡張子はいらない
                    "objfile": objfolder + "/tank.obj"
                    }

        object = {"name": "object",
                  "type": "RigidBody",
                  "output_vtu_file_name": "floatingbody",  # 拡張子はいらない
                  "output_pvd_file_name": "floatingbody",  # 拡張子はいらない
                  "velocity": "floating",
                  "objfile": objfolder + "/object.obj"}

        input_files = [wavetank, water]

        setting = {"RK_order": 1,
                   "max_dt": 0.001,
                   "CSML": 3.05,
                   "end_time_step": 10000,
                   "end_time": 10,
                   "initial_surface_z_position": 0.1,
                   # "particle_spacing": 0.00625,
                   "particle_spacing": 0.2/25,
                   "input_files": [x["name"]+".json" for x in input_files]}
    case "Lobovsky2014":

        objfolder = program_home + "/cpp/obj/2022Arai/Lobovsky2014"

        water = {"name": "water",
                 "type": "Fluid",
                 "output_vtu_file_name": "water",  # 拡張子はいらない
                 "output_pvd_file_name": "water",  # 拡張子はいらない
                 "objfile": objfolder + "/water.obj"}

        wavetank = {"name": "wavetank",
                    "type": "RigidBody",
                    # "ignore": wavetank_ignore,
                    "output_vtu_file_name": "wavetank",  # 拡張子はいらない
                    "output_pvd_file_name": "wavetank",  # 拡張子はいらない
                    "objfile": objfolder + "/tank5.obj"}

        sensor1 = {"name": "sensor1",
                   "type": "probe",
                   "location": [1610 / 1000., 150 / 2. / 1000., 3 / 1000.]}

        sensor2 = {"name": "sensor2",
                   "type": "probe",
                   "location": [1610 / 1000., 150 / 2. / 1000., 15 / 1000.]}

        sensor2L = {"name": "sensor2L",
                    "type": "probe",
                    "location": [1610 / 1000., (150 - 37.5) / 2. / 1000., 15 / 1000.]}

        sensor3 = {"name": "sensor3",
                   "type": "probe",
                   "location": [1610 / 1000., 150 / 2. / 1000., 30 / 1000.]}

        sensor4 = {"name": "sensor4",
                   "type": "probe",
                   "location": [1610 / 1000., 150 / 2. / 1000., 80 / 1000.]}

        input_files = [wavetank, water,
                       sensor1, sensor2, sensor2L, sensor3, sensor4]

        setting = {"RK_order": 1,
                   "max_dt": 0.001,
                   "end_time_step": 20000,
                   "end_time": 1,
                   "CSML": 3.05,
                   "initial_surface_z_position": 0.6,
                   "particle_spacing": 0.015,
                   "input_files": [x["name"]+".json" for x in input_files]}
    case "Kamra2019":

        id = "_square"

        objfolder = program_home + "/cpp/obj/2022Arai/Kamra2019"

        water = {"name": "water",
                 "type": "Fluid",
                 "output_vtu_file_name": "water",  # 拡張子はいらない
                 "output_pvd_file_name": "water",  # 拡張子はいらない
                 "objfile": objfolder + "/water.obj"}

        wavetank = {"name": "wavetank",
                    "type": "RigidBody",
                    # "ignore": wavetank_ignore,
                    "output_vtu_file_name": "wavetank",  # 拡張子はいらない
                    "output_pvd_file_name": "wavetank",  # 拡張子はいらない
                    "objfile": objfolder + "/tank"+id+"_cylinder.obj"}

        sensor1 = {"name": "sensor1",
                   "type": "probe",
                   "location": [0.6, 0.1, 0.011]}

        sensor2 = {"name": "sensor2",
                   "type": "probe",
                   "location": [0.8, 0.1, 0.004]}

        input_files = [wavetank, water,  sensor1, sensor2]

        setting = {"RK_order": 2,
                   "max_dt": 0.0001,
                   "end_time_step": 50000,
                   "end_time": 0.5,
                   "CSML": 2.7,
                   "initial_surface_z_position": 0.2,
                   "particle_spacing": 0.0125,
                   "input_files": [x["name"]+".json" for x in input_files]}

# ---------------------------------------------------------------------------- #

id = id + "_PS" + str(setting["particle_spacing"]).replace(".", "d") \
    + "_CSML" + str(setting["CSML"]).replace(".", "d")\
    + "_RK" + str(setting["RK_order"])
input_directory = "./input_files/" + SimulationCase + id
output_directory = home + "/SPH/" + SimulationCase + id
os.makedirs(input_directory, exist_ok=True)
os.makedirs(output_directory, exist_ok=True)
setting["output_directory"] = output_directory

# ---------------------------------------------------------------------------- #

red = '\033[91m'
green = '\033[92m'
magenta = '\033[95m'
coloroff = '\033[0m'
# @ -------------------------------------------------------- #
# @           その他，water.json,tank.json などを出力           #
# @ -------------------------------------------------------- #
for INPUTS in input_files:
    print('------------------------------------')
    for key, value in INPUTS.items():
        print(f'{key: <{20}}', '\t', green, value, coloroff)
    print('------------------------------------')
    f = open(input_directory+"/"+INPUTS["name"]+".json", 'w')
    json.dump(INPUTS, f, ensure_ascii=True, indent=4)
    f.close()

# @ -------------------------------------------------------- #
# @                  setting.json を出力                      #
# @ -------------------------------------------------------- #
print('------------------------------------')
for key, value in setting.items():
    print(f'{key: <{20}}', '\t', red, value, coloroff)
print('------------------------------------')
f = open(input_directory+"/setting.json", 'w')
json.dump(setting, f, ensure_ascii=True, indent=4)
f.close()

print("The directory for input files :", magenta, input_directory, coloroff)
