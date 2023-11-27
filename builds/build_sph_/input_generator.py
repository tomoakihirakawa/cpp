
from os.path import expanduser
import os
import math
import json
from math import pi
import platform

import sys
sys.path.append("..")
from IOgenerator import generate_input_files

current_directory = os.path.dirname(os.path.abspath(__file__))
program_home = os.path.join(current_directory, '../../../../code')

rho = 1000.
g = 9.81

# sh clean
# cmake -DCMAKE_BUILD_TYPE=Release ../
# make


# ---------------------------------------------------------------------------- #

def IO_dir(id):
    home = expanduser("~")
    input_directory = "./input_files/" + id
    os.makedirs(input_directory, exist_ok=True)
    output_directory = home + "/SPH/" + id
    os.makedirs(output_directory, exist_ok=True)
    return input_directory, output_directory

SimulationCase = "Lobovsky2013"
id = ""
match SimulationCase:
    case "static_pressure":

        objfolder = program_home + "/cpp/obj/SPH/open_closed_tank"

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

        particle_spacing = 0.2/10

        setting = {"RK_order": 1,
                   "max_dt": particle_spacing/10,
                   "CSML": 2.8,
                   "end_time_step": 1000000,
                   "end_time": 10,
                   "initial_surface_z_position": 0.1,
                   "particle_spacing": particle_spacing}

        id = SimulationCase + "_PS" + str(setting["particle_spacing"]).replace(".", "d") \
                            + "_CSML" + str(setting["CSML"]).replace(".", "d")\
                            + "_RK" + str(setting["RK_order"])

        generate_input_files(input_files, setting,IO_dir, id)

    case "Lobovsky2013":

        objfolder = program_home + "/cpp/obj/Lobovsky2013/original"

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
                   "max_dt": 0.0002,
                   "end_time_step": 20000,
                   "end_time": 1,
                   "CSML": 2.5,
                   "initial_surface_z_position": 0.6,
                   "particle_spacing": 0.02}

        id = SimulationCase + "_PS" + str(setting["particle_spacing"]).replace(".", "d") \
                            + "_CSML" + str(setting["CSML"]).replace(".", "d")\
                            + "_RK" + str(setting["RK_order"])

        generate_input_files(input_files, setting, IO_dir, id)

    case "Kamra2019":

        id = "_square"

        objfolder = program_home + "/cpp/obj/Kamra2019"

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
                   "end_time": 4.,
                   "CSML": 2.4,
                   "initial_surface_z_position": 0.2,
                   "particle_spacing": 0.015}

        id = SimulationCase + id \
                            + "_PS" + str(setting["particle_spacing"]).replace(".", "d") \
                            + "_CSML" + str(setting["CSML"]).replace(".", "d")\
                            + "_RK" + str(setting["RK_order"])

        generate_input_files(input_files, setting, IO_dir, id)
