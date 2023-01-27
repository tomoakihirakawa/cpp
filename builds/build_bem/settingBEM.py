from os.path import expanduser
import os
import math
import json
from math import pi
import platform

home = expanduser("~")

if platform.system() == "Linux":
    program_home = home + "/research/"
else:
    program_home = home + "/Dropbox/markdown/"

rho = 1000.
g = 9.81

'''
プログラムを回す際に面倒な事は，入力ファイルの設定方法．
入力ファイルの作り方をドキュメントで示されても，具体的な例がないとわかりにくい．
例があっても，例と違う場合どうすればいいかなど，わからないことは多い．
このように，入力ファイルを生成するプログラムを作っておけば，その面倒をだいぶ解消できる．
'''

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

SimulationCase = "three_floatingbodies"
match SimulationCase:
    case "2022Tsukada_flotingbody_without_moonpool_A0d75T7d0":

        input_directory = "./inputs_" + SimulationCase
        os.makedirs(input_directory, exist_ok=True)
        output_directory = home + "/BEM/2022Tsukada_flotingbody_without_moonpool_A0d75T7d0"
        os.makedirs(output_directory, exist_ok=True)

        water = {"name": "water", "type": "Fluid"}

        tank = {"name": "tank", "type": "RigidBody", "isFixed": True}

        start = 0.
        a = 2.
        T = 7.
        h = 80
        z_surface = 80

        wavemaker = {"name": "wavemaker",
                     "type": "SoftBody",
                     "isFixed": True,
                     "velocity": ["linear_traveling_wave", start, a, T, h, z_surface]}

        floatingbody = {"name": "floatingbody",
                        "type": "RigidBody",
                        "velocity": "floating"}

        A = 2450.00
        floatingbody["mass"] = m = (1000.*g*7.5*A)/g
        floatingbody["COM"] = [100., 50., 75.]
        floatingbody["radius_of_gyration"] = [20., 20., 20.]
        floatingbody["MOI"] = [m*math.pow(floatingbody["radius_of_gyration"][0], 2),
                               m*math.pow(floatingbody["radius_of_gyration"][1], 2),
                               m*math.pow(floatingbody["radius_of_gyration"][2], 2)]

        objfolder = program_home + "/cpp/obj/tsukada2022_re"
        water["objfile"] = objfolder + "/water100_mod.obj"
        wavemaker["objfile"] = objfolder + "/wavemaker200.obj"
        tank["objfile"] = objfolder + "/wavetank10.obj"
        floatingbody["objfile"] = objfolder+"/floatingbody100.obj"

        inputfiles = [tank, wavemaker, water, floatingbody]

        setting = {"max_dt": 0.2,
                   "end_time_step": 10000,
                   "end_time": 10000,
                   "output_directory": output_directory,
                   "input_files": [x["name"]+".json" for x in inputfiles]}
    case "Tonegawa20221216":

        input_directory = "./inputs3"
        os.makedirs(input_directory, exist_ok=True)
        output_directory = home + "/BEM/2022Tonegawa3"
        os.makedirs(output_directory, exist_ok=True)

        water = {"name": "water",
                 "type": "Fluid"
                 }

        tank = {"name": "tank",
                "type": "RigidBody",
                "isFixed": True}

        start = 0.
        a = 1.
        T = 7.
        h = 150
        z_surface = 150

        wavemaker = {"name": "wavemaker",
                     "type": "SoftBody",
                     "isFixed": True,
                     "velocity": ["linear_traveling_wave", start, a, T, h, z_surface]}

        floatingbody = {"name": "floatingbody",
                        "type": "RigidBody",
                        "velocity": "floating"}

        A = 162.122
        floatingbody["mass"] = m = (1000.*g*78*A)/g
        floatingbody["COM"] = [75., 50., 150-78./2]
        floatingbody["radius_of_gyration"] = [20., 20., 20.]
        floatingbody["MOI"] = [m*math.pow(floatingbody["radius_of_gyration"][0], 2),
                               m*math.pow(floatingbody["radius_of_gyration"][1], 2),
                               m*math.pow(floatingbody["radius_of_gyration"][2], 2)]

        objfolder = home + "/Dropbox/markdown/cpp/obj/2022Tonegawa/test20221218"
        water["objfile"] = objfolder + "/water150_mod.obj"
        wavemaker["objfile"] = objfolder + "/wavemaker100.obj"
        tank["objfile"] = objfolder + "/tank10.obj"
        floatingbody["objfile"] = objfolder+"/floatingbody50.obj"

        inputfiles = [tank, wavemaker, water, floatingbody]

        setting = {"max_dt": 0.2,
                   "end_time_step": 10000,
                   "end_time": 10000,
                   "output_directory": output_directory,
                   "input_files": [x["name"]+".json" for x in inputfiles]}
    case "two_floatingbodies":

        start = 0.
        a = 1.5
        T = 7.  # 5-8
        h = 150
        z_surface = 150

        id = "_a" + str(a).replace(".", "d")\
            + "_T" + str(T).replace(".", "d")\
            + "_h" + str(h).replace(".", "d")

        input_directory = "./inputs_" + SimulationCase + id
        os.makedirs(input_directory, exist_ok=True)
        output_directory = home + "/BEM/" + SimulationCase + id
        os.makedirs(output_directory, exist_ok=True)

        water = {"name": "water",
                 "type": "Fluid"}

        tank = {"name": "tank",
                "type": "RigidBody",
                "isFixed": True}

        wavemaker = {"name": "wavemaker",
                     "type": "SoftBody",
                     "isFixed": True,
                     "velocity": ["linear_traveling_wave", start, a, T, h, z_surface]}

        floatingbody_a = {"name": "floatingbody_a",
                          "COM": [300., 150., 150-78./2]}

        floatingbody_b = {"name": "floatingbody_b",
                          "COM": [600., 150., 150-78./2]}

        floatingbody_c = {"name": "floatingbody_c",
                          "COM": [900., 150., 150-78./2]}

        A = 170.779
        for x in [floatingbody_a, floatingbody_b, floatingbody_c]:
            x["type"] = "RigidBody"
            x["velocity"] = "floating"
            x["mass"] = m = (1000.*g*78*A)/g
            x["radius_of_gyration"] = [20., 20., 20.]
            x["MOI"] = [m*math.pow(x["radius_of_gyration"][0], 2),
                        m*math.pow(x["radius_of_gyration"][1], 2),
                        m*math.pow(x["radius_of_gyration"][2], 2)]

        # ------------------------------------------------------------------------#

        objfolder = program_home + "/cpp/obj/2022Tonegawa/two_floatingbodies/"
        water["objfile"] = objfolder + "/water300_two_mod2.obj"
        wavemaker["objfile"] = objfolder + "/wavemaker100_two.obj"
        tank["objfile"] = objfolder + "/tank10_two.obj"
        floatingbody_a["objfile"] = objfolder + "floatingbody_a50.obj"
        floatingbody_b["objfile"] = objfolder + "floatingbody_b50_two.obj"
        floatingbody_c["objfile"] = objfolder + "floatingbody_c50_two.obj"

        inputfiles = [tank, wavemaker, water,
                      floatingbody_b, floatingbody_c]

        setting = {"max_dt": 0.15,
                   "end_time_step": 12000,
                   "end_time": 120,
                   "output_directory": output_directory,
                   "input_files": [x["name"]+".json" for x in inputfiles]}
    case "three_floatingbodies":

        start = 0.
        a = 1.5
        T = 5.0  # 5-8
        h = 150
        z_surface = 150

        id = "_a" + str(a).replace(".", "d")\
            + "_T" + str(T).replace(".", "d")\
            + "_h" + str(h).replace(".", "d")

        input_directory = "./inputs_" + SimulationCase + id
        os.makedirs(input_directory, exist_ok=True)
        output_directory = home + "/BEM/" + SimulationCase + id
        os.makedirs(output_directory, exist_ok=True)

        water = {"name": "water",
                 "type": "Fluid"}

        tank = {"name": "tank",
                "type": "RigidBody",
                "isFixed": True}

        wavemaker = {"name": "wavemaker",
                     "type": "SoftBody",
                     "isFixed": True,
                     "velocity": ["linear_traveling_wave", start, a, T, h, z_surface]}

        floatingbody_a = {"name": "floatingbody_a",
                          "COM": [300., 150., 150-78./2]}

        floatingbody_b = {"name": "floatingbody_b",
                          "COM": [600., 150., 150-78./2]}

        floatingbody_c = {"name": "floatingbody_c",
                          "COM": [900., 150., 150-78./2]}

        A = 170.779
        for x in [floatingbody_a, floatingbody_b, floatingbody_c]:
            x["type"] = "RigidBody"
            x["velocity"] = "floating"
            x["mass"] = m = (1000.*g*78*A)/g
            x["radius_of_gyration"] = [20., 20., 20.]
            x["MOI"] = [m*math.pow(x["radius_of_gyration"][0], 2),
                        m*math.pow(x["radius_of_gyration"][1], 2),
                        m*math.pow(x["radius_of_gyration"][2], 2)]

        # ------------------------------------------------------------------------#

        objfolder = program_home + "/cpp/obj/2022Tonegawa/float_case1/"
        water["objfile"] = objfolder + "/water300_mod.obj"
        wavemaker["objfile"] = objfolder + "/wavemaker100.obj"
        tank["objfile"] = objfolder + "/tank10.obj"
        floatingbody_a["objfile"] = objfolder + "floatingbody_a50.obj"
        floatingbody_b["objfile"] = objfolder + "floatingbody_b50.obj"
        floatingbody_c["objfile"] = objfolder + "floatingbody_c50.obj"

        inputfiles = [tank, wavemaker, water,
                      floatingbody_a, floatingbody_b, floatingbody_c]

        setting = {"max_dt": 0.15,
                   "end_time_step": 12000,
                   "end_time": 120,
                   "output_directory": output_directory,
                   "input_files": [x["name"]+".json" for x in inputfiles]}
    case "Retzler2000simple":

        input_directory = "./inputs_" + SimulationCase
        os.makedirs(input_directory, exist_ok=True)
        output_directory = home + "/BEM/" + SimulationCase
        os.makedirs(output_directory, exist_ok=True)

        water = {"name": "water",
                 "type": "Fluid",
                 }

        tank = {"name": "tank",
                "type": "RigidBody",
                "isFixed": True}

        wavemaker = {"name": "wavemaker",
                     "type": "RigidBody",
                     "isFixed": True,
                     "velocity": ["Retzler2000", 0.1]}

        objfolder = program_home + "/cpp/obj/chaplin2000simple/"
        water["objfile"] = objfolder+"/water_refined.obj"
        wavemaker["objfile"] = objfolder+"/cylinder200.obj"
        tank["objfile"] = objfolder+"/tank.obj"
        wavetank_ignore = False
        inputfiles = [tank, wavemaker, water]
        setting = {"max_dt": 0.002,
                   "end_time_step": 1000,
                   "end_time": 0.4,
                   "output_directory": output_directory,
                   "input_files": [x["name"]+".json" for x in inputfiles]}
    case "Retzler2000":

        input_directory = "./inputs_Retzler2000simple"
        os.makedirs(input_directory, exist_ok=True)
        output_directory = home + "/BEM/Retzler2000simple"
        os.makedirs(output_directory, exist_ok=True)

        water = {"name": "water",
                 "type": "Fluid",
                 }

        tank = {"name": "tank",
                "type": "RigidBody",
                "isFixed": True}

        wavemaker = {"name": "wavemaker",
                     "type": "RigidBody",
                     "isFixed": True,
                     "velocity": ["Retzler2000", 0.1]}

        objfolder = home + "/Dropbox/markdown/cpp/obj/chaplin2000/"
        water["objfile"] = objfolder+"/water_remeshed2.obj"
        wavemaker["objfile"] = objfolder+"/cylinder200.obj"
        tank["objfile"] = objfolder+"/tank.obj"
        wavetank_ignore = False
        inputfiles = [tank, wavemaker, water]
        setting = {"max_dt": 0.002,
                   "end_time_step": 100,
                   "end_time": 0.4,
                   "output_directory": output_directory,
                   "input_files": [x["name"]+".json" for x in inputfiles]}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
white = '\033[90m'
red = '\033[91m'
blue = '\033[96m'
green = '\033[92m'
magenta = '\033[95m'
coloroff = '\033[0m'
#@ -------------------------------------------------------- #
#@                  setting.json を出力                      #
#@ -------------------------------------------------------- #
print('------------------------------------')
for key, value in setting.items():
    print(f'{key: <{20}}', '\t', green, value, coloroff)
print('------------------------------------')
f = open(input_directory+"/setting.json", 'w')
json.dump(setting, f, ensure_ascii=True, indent=4)
f.close()

#@ -------------------------------------------------------- #
#@           その他，water.json,tank.json などを出力           #
#@ -------------------------------------------------------- #
print("The directory for input files :", magenta, input_directory, coloroff)
for INPUTS in inputfiles:
    print('------------------------------------')
    for key, value in INPUTS.items():
        if value == "floating":
            print(f'{key: <{20}}', '\t', green, value, coloroff)
        elif value == "RigidBody":
            print(f'{key: <{20}}', '\t', red, value, coloroff)
        elif value == "Fluid":
            print(f'{key: <{20}}', '\t', blue, value, coloroff)
        else:
            print(f'{key: <{20}}', '\t', white, value, coloroff)
    print('------------------------------------')
    f = open(input_directory+"/"+INPUTS["name"]+".json", 'w')
    json.dump(INPUTS, f, ensure_ascii=True, indent=4)
    f.close()
