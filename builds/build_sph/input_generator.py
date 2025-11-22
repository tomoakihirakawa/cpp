

import copy
import platform
import json
import math
import os
import sys
import argparse

from math import pi
from os.path import expanduser

import sys
sys.path.append("..")
input_dir = "./input_files/"
sys_home_dir = expanduser("~")
current_dir = os.path.dirname(os.path.abspath(__file__))
code_home_dir = os.path.join(current_dir, '../../../')

from IOgenerator import generate_input_files, read_obj, Norm3d, Differece3d, Add3d

rho = 1000.
g = 9.81

# ---------------------------------------------------------------------------- #

parser = argparse.ArgumentParser(description='コマンドライン引数を受け取るためのパーサー')
parser.add_argument('case', type=str, help='シミュレーションケース')
parser.add_argument('-m', '--mesh', type=str, help='メッシュの名前')
parser.add_argument('-wavemaker', type=str, help='波を作る方法')
parser.add_argument('-dt', '--max_dt', type=float, required=True, help='時間刻み幅')
parser.add_argument('-H', '--wave_height', type=float, help='波の高さ')
parser.add_argument('-o', '--outputdir', type=str, help='出力ディレクトリ')
parser.add_argument('-T', '--wave_period', type=float, help='波の周期')
parser.add_argument('-L', '--wave_length', type=float, help='波の波長')
parser.add_argument('-theta', type=float, help='波の進行方向')
parser.add_argument('-s', '--suffix', type=str, help='接尾語')

args = parser.parse_args()

# SimulationCase = ""
meshname = ""
wavemaker_type = ""
suffix = ""
default_output_dir = None
L = None
padding = 20

SimulationCase = args.case
print(f"{'SimulationCase:':>{padding}} {SimulationCase}")

if args.mesh is not None:
    meshname = args.mesh
    print(f"meshname: {meshname}")
if args.wavemaker is not None:
    wavemaker_type = args.wavemaker
    print(f"wavemaker_type: {wavemaker_type}")
if args.max_dt:
    dt = args.max_dt
    print(f"{'max_dt:':>{padding}} {dt}")
if args.wave_period is not None:
    T = args.wave_period
    print(f"{'T:':>{padding}} {T}")
if args.wave_length is not None:
    L = args.wave_length
    print(f"{'L:':>{padding}} {L}")
if args.suffix is not None:
    suffix = args.suffix
    print(f"{'suffix:':>{padding}} {suffix}")
if args.wave_height is not None:
    H = args.wave_height
    print(f"{'H:':>{padding}} {H}")
if args.outputdir is not None:
    default_output_dir = args.outputdir
    print(f"{'outputdir:':>{padding}} {default_output_dir}")
if args.theta is not None:
    theta = args.theta
    print(f"{'theta:':>{padding}} {theta}")
else:
    theta = None

# ---------------------------------------------------------------------------- #

def IO_dir(id):
    home = expanduser("~")
    input_directory = "./input_files/" + id
    # home = "/Volumes/home"
    os.makedirs(input_directory, exist_ok=True)    
    output_directory = home + "/SPH/" + id
    os.makedirs(output_directory, exist_ok=True)
    return input_directory, output_directory

SimulationCase = "static_pressure"
id = ""

if len(sys.argv) > 1:
    # If a value is passed, use it
    SimulationCase = sys.argv[1]
    print(f"Value passed from command line: {SimulationCase}")
else:
    # Default action if no value is passed
    print("No value passed. Executing default action.")

if "static_pressure" in SimulationCase:

    objfolder = code_home_dir + "/cpp/obj/SPH/open_closed_tank"

    water = {"name": "water",
                "type": "Fluid",
                "objfile": objfolder + "/water.obj"}

    wavetank = {"name": "wavetank",
                "type": "RigidBody",
                "objfile": objfolder + "/open_tank.obj"}

    object = {"name": "object",
                "type": "RigidBody",
                "velocity": "floating",
                "objfile": objfolder + "/object.obj"}

    input_files = [wavetank, water]

    particle_spacing = 0.01

    setting = {"RK_order": 1,
                "max_dt": particle_spacing/5,
                "CSML": 2.4,
                "end_time_step": 1000000,
                "end_time": 100,
                "initial_surface_z_position": 0.1,
                "particle_spacing": particle_spacing}

    id = SimulationCase + "_PS" + str(setting["particle_spacing"]).replace(".", "d") 
    # \                   + "_WITHOUT_CORRECTION"

    generate_input_files(input_files, setting,IO_dir, id)

elif "Lobovsky2013_small" in SimulationCase:

    objfolder = code_home_dir + "/cpp/obj/Lobovsky2013_small"

    H = 0.3
    
    water = {"name": "water",
                "type": "Fluid",
                "objfile": objfolder + "/water"+"H"+str(H).replace(".", "d")+".obj"}

    gate = {"name": "gate",
            "type": "RigidBody",
            "inactivate" : [0.001, 1000.],
            "objfile": objfolder + "/gate.obj"}

    wavetank = {"name": "wavetank",
                "type": "RigidBody",
                "objfile": objfolder + "/tank.obj"}

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

    input_files = [wavetank, water, gate, sensor1, sensor2, sensor2L, sensor3, sensor4]

    setting = {"RK_order": 1,
                "max_dt": 0.0004,
                "end_time_step": 100000*5,
                "end_time": 5,
                "CSML": 2.8,
                "initial_surface_z_position": 0.6,
                "particle_spacing": 0.015}

    id = SimulationCase + "_PS" + str(setting["particle_spacing"]).replace(".", "d") \
                        + "_CSML" + str(setting["CSML"]).replace(".", "d")\
                        + "_RK" + str(setting["RK_order"])

    generate_input_files(input_files, setting, IO_dir, id)

elif "Lobovsky2013" in SimulationCase:

    objfolder = code_home_dir + "/cpp/obj/Lobovsky2013/original"

    H = 0.3
    
    water = {"name": "water",
                "type": "Fluid",
                "objfile": objfolder + "/water"+"H"+str(H).replace(".", "d")+".obj"}

    gate = {"name": "gate",
            "type": "RigidBody",
            "inactivate" : [0.001, 1000.],
            "objfile": objfolder + "/gate.obj"}

    wavetank = {"name": "wavetank",
                "type": "RigidBody",
                "objfile": objfolder + "/tank.obj"}

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

    input_files = [wavetank, water, gate, sensor1, sensor2, sensor2L, sensor3, sensor4]

    setting = {"RK_order": 1,
                "max_dt": dt,
                "end_time_step": 3000*5,
                "end_time": 2,
                "CSML": 2.5,
                "initial_surface_z_position": 0.3,
                "particle_spacing": 0.018}

    id = SimulationCase + "_PS" + str(setting["particle_spacing"]).replace(".", "d") \
                        + "_CSML" + str(setting["CSML"]).replace(".", "d")\
                        + "_RK" + str(setting["RK_order"])

    generate_input_files(input_files, setting, IO_dir, id)

elif "Kamra2019" in SimulationCase:

    id = "_square"

    objfolder = code_home_dir + "/cpp/obj/Kamra2019"

    water = {"name": "water",
                "type": "Fluid",
                "objfile": objfolder + "/water.obj"}

    wall = {"name": "wall",
            "type": "RigidBody",
            "inactivate" : [0.0001, 1000.],
            "objfile": objfolder + "/wall.obj"}

    wavetank = {"name": "tank",
                "type": "RigidBody",
                "objfile": objfolder + "/tank"+id+".obj"}

    sensor1 = {"name": "sensor1",
                "type": "probe",
                "location": [0.6, 0.1, 0.011]}

    sensor2 = {"name": "sensor2",
                "type": "probe",
                "location": [0.8, 0.1, 0.004]}

    # input_files = [wavetank, water, wall, sensor1, sensor2]

    input_files = [wavetank, water, sensor1, sensor2]

    setting = {"RK_order": 1,  # \label{SPH:RK_order}
                "max_dt": 0.0005,
                "end_time_step": 100000,
                "end_time": 10.,
                "CSML": 2.5,
                "initial_surface_z_position": 0.2,
                "particle_spacing": 0.0125}

    id = SimulationCase + id \
                        + "_PS" + str(setting["particle_spacing"]).replace(".", "d") \
                        + "_CSML" + str(setting["CSML"]).replace(".", "d")\
                        + "_RK" + str(setting["RK_order"])

    generate_input_files(input_files, setting, IO_dir, id)
