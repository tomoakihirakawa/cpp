from math import pi
import json
import math
from os.path import expanduser
home = expanduser("~")

# -------------------------------------------------------- #
# data = {
#     "output_file_name": "water.pvd",
#     "output_directory": "BEM",
#     "radius": 0.05,
# }

# f = open("./settingBEM.json", 'w')
# json.dump(data, f, ensure_ascii=True, indent=4)
# f.close()

#! -------------------------------------------------------- #

# output_directory = home+"/BEM/BEM_dt0d01_sloshing_Xdir_H0d10_L0d25_A0d01_ww1_1d0_using_Dombre2019_EMT_K0d0"
# output_directory = home+"/BEM/test"
settingBEM = {
    "max_dt": 0.001,
    "mesh": 90,
    "stop_remesh_time": 150.,
    "force_remesh_time": 0.1,
    "preparation_time": 0.0,
    "preparation_max_dt": 0.001,
    "grid_refinement":  0
}

#! -------------------------------------------------------- #
FloatingWind1 = False
sloshing_tank = False
sloshing_tank_torus = False
sloshing_tank_sphere = False
Li2002 = False
Chaplin2000 = False
Goring1979 = False
Chaplin2000simple = False
Retzler2000simple = False
Retzler2000 = False
#! -------------------------------------------------------- #
Goring1979 = True
#! -------------------------------------------------------- #
if FloatingWind1:
    settingBEM["max_dt"] = 0.02
    objfolder = "FloatingWind1"
    water = "/water60.obj"
    floatingbody = "floatingbody.obj"
    wavemaker = "tank.obj"
    wavetank_ignore = False
    settingBEM["output_directory"] = home + "/BEM/FloatingWind1"
    settingBEM["inputfiles"] = ["water.json",
                                "floatingbody.json", "wavemaker.json"]
elif sloshing_tank_torus:
    # water = "/water"+str(settingBEM["mesh"])+".obj"
    water = "/water60.obj"
    objfolder = "sloshing_tank_torus"
    wavemaker = "torus4.obj"
    wavetank = "tank3.obj"
    wavetank_ignore = False
    settingBEM["output_directory"] = home + "/BEM/tourus"
    settingBEM["inputfiles"] = ["water.json",
                                "floatingbody.json", "wavetank.json"]
elif sloshing_tank_sphere:
    # water = "/water"+str(settingBEM["mesh"])+".obj"
    water = "/water20.obj"
    objfolder = "sloshing_tank_sphere"
    wavemaker = "sphere4.obj"
    #
    wavetank = "tank3.obj"
    wavetank_ignore = False
    settingBEM["output_directory"] = home + "/BEM/sphere"
    settingBEM["inputfiles"] = ["water.json",
                                "wavemaker.json", "wavetank.json"]
elif sloshing_tank:
    # water = "/water"+str(settingBEM["mesh"])+".obj"
    water = "/water20.obj"
    objfolder = "sloshing_tank"
    wavemaker = "tank3.obj"
    wavetank = "tank3.obj"
    wavetank_ignore = True
    settingBEM["output_directory"] = home + "/BEM/sloshing_tank2"
    settingBEM["inputfiles"] = ["water.json",
                                "wavemaker.json", "wavetank.json"]
elif Retzler2000simple:
    objfolder = "chaplin2000simple"
    water = "/water.obj"
    wavemaker = "cylinder100.obj"
    wavetank = "tank.obj"
    wavetank_ignore = False
    settingBEM["inputfiles"] = ["water.json",
                                "wavemaker.json", "wavetank.json"]
    settingBEM["max_dt"] = 0.002
    settingBEM["output_directory"] = home + "/BEM/Retzler2000_" + \
        str(settingBEM["max_dt"]).replace(
            '.', 'd') + "AreaWeightedVelocity_1_Linear_simple"
    # output_directory = home + "/BEM/Retzler2000_0d005_linear"
    # output_directory = home + "/BEM/test"
elif Goring1979:
    mesh = "150_refined"
    water = "/water"+mesh+".obj"
    objfolder = "/Goring1979"
    wavemaker = "wavemaker.obj"
    wavetank = "wavetank.obj"
    wavetank_ignore = False
    settingBEM["inputfiles"] = ["water.json",
                                "wavemaker.json", "wavetank.json"]
    settingBEM["max_dt"] = 0.04
    settingBEM["output_directory"] = home + "/BEM/Goring1979_" + \
        str(settingBEM["max_dt"]).replace('.', 'd') + \
        "_Linear" + mesh + "_1times"
elif Li2002:
    water = "/water60_refined.obj"
    objfolder = "/Li2002"
    wavemaker = "wavemaker.obj"
    wavetank = "tank.obj"
    wavetank_ignore = False
    settingBEM["inputfiles"] = ["water.json",
                                "wavemaker.json", "wavetank.json"]
    settingBEM["max_dt"] = 0.01
    settingBEM["output_directory"] = home + "/BEM/Li2002_" + \
        str(settingBEM["max_dt"]).replace('.', 'd') + \
        "_Hh0d2_AreaWeightedVelocity_1_Linear1"
    # output_directory = home + "/BEM/Retzler2000_0d005_linear"
    # output_directory = home + "/BEM/test"
elif Retzler2000:
    water = "/water_remeshed2.obj"
    objfolder = "chaplin2000"
    wavemaker = "cylinder100.obj"
    wavetank = "tank3.obj"
    wavetank_ignore = False
    settingBEM["inputfiles"] = ["water.json",
                                "wavemaker.json", "wavetank.json"]
    settingBEM["max_dt"] = 0.002
    settingBEM["output_directory"] = home + "/BEM/Retzler2000_" + \
        str(settingBEM["max_dt"]).replace('.', 'd') + \
        "_Linear_10times"
    # output_directory = home + "/BEM/Retzler2000_0d005_linear"
    # output_directory = home + "/BEM/test"
elif Chaplin2000:
    # water = "/water30.obj"
    # water = "/water115_test.obj"
    water = "/water_remeshed2.obj"
    # water = "/water80_refined_120.obj"
    # water = "/water80_test.obj"
    # water = "/water100_test.obj"
    # water = "/water100_refined_80.obj"
    # water = "/water30_refined_199.obj"
    # water = "/water30_refined_91.obj"
    objfolder = "chaplin2000"
    wavemaker = "cylinder100.obj"
    wavetank = "tank3.obj"
    wavetank_ignore = False
    settingBEM["inputfiles"] = ["water.json",
                                "wavemaker.json", "wavetank.json"]
    settingBEM["max_dt"] = 0.01
    settingBEM["output_directory"] = home + "/BEM/Chaplin2000_" + \
        str(settingBEM["max_dt"]).replace('.', 'd') + \
        "_Linear_10times"
    # output_directory = home + "/BEM/Retzler2000_0d005_linear"
    # output_directory = home + "/BEM/test"
elif Chaplin2000simple:
    objfolder = "chaplin2000simple"
    water = "/water.obj"
    wavemaker = "cylinder100.obj"
    wavetank = "tank.obj"
    wavetank_ignore = False
    settingBEM["inputfiles"] = ["water.json",
                                "wavemaker.json",
                                "wavetank.json"]
    settingBEM["max_dt"] = 0.005
    settingBEM["output_directory"] = home + "/BEM/test__" + \
        str(settingBEM["max_dt"]).replace(
            '.', 'd') + "_linear3_"
elif refinement:
    water = "/water_remeshed2.obj"
    objfolder = "chaplin2000"
    wavemaker = "cylinder100.obj"
    wavetank = "tank3.obj"
    wavetank_ignore = False
    settingBEM["output_directory"] = home + "/BEM/refinement"
    settingBEM["grid_refinement"] = 200
    settingBEM["inputfiles"] = ["water.json",
                                "wavemaker.json", "wavetank.json"]

#! -------------------------------------------------------- #

print('------------------------------------')
for key, value in settingBEM.items():
    print(key, '\t', value)
print('------------------------------------')
f = open("./setting.json", 'w')
json.dump(settingBEM, f, ensure_ascii=True, indent=4)
f.close()

#! -------------------------------------------------------- #

data = {
    "name": "water",
    "objfile": "../../obj/"+objfolder+water,
    "type": "Fluid",
    "output_vtu_file_name": "water",  # 拡張子はいらない
    "output_pvd_file_name": "water",  # 拡張子はいらない
}

print('------------------------------------')
for key, value in data.items():
    print(key, '\t', value)
print('------------------------------------')
f = open("./water.json", 'w')
json.dump(data, f, ensure_ascii=True, indent=4)
f.close()

#! -------------------------------------------------------- #
try:
    data = {
        "name": "wavemaker",
        "objfile": "../../obj/"+objfolder+"/"+wavemaker,
        "type": "RigidBody",
        "output_vtu_file_name": "wavemaker",  # 拡張子はいらない
        "output_pvd_file_name": "wavemaker",  # 拡張子はいらない
        "velocity": ["Goring1979", 4.],
    }

    print('------------------------------------')
    for key, value in data.items():
        print(key, '\t', value)
    print('------------------------------------')
    f = open("./wavemaker.json", 'w')
    json.dump(data, f, ensure_ascii=True, indent=4)
    f.close()
except:
    print("no wavemaker")
#! -------------------------------------------------------- #
try:
    data = {
        "name": "wavetank",
        "objfile": "../../obj/"+objfolder+"/"+wavetank,
        "type": "RigidBody",
        "ignore": wavetank_ignore,
        "output_vtu_file_name": "wavetank",  # 拡張子はいらない
        "output_pvd_file_name": "wavetank",  # 拡張子はいらない
    }

    print('------------------------------------')
    for key, value in data.items():
        print(key, '\t', value)
    print('------------------------------------')
    f = open("./wavetank.json", 'w')
    json.dump(data, f, ensure_ascii=True, indent=4)
    f.close()
except:
    print("no wavetank")
#! -------------------------------------------------------- #
try:
    data = {
        "name": "floatingbody",
        "objfile": "../../obj/"+objfolder+"/"+floatingbody,
        "type": "RigidBody",
        "output_vtu_file_name": "floatingbody",  # 拡張子はいらない
        "output_pvd_file_name": "floatingbody",  # 拡張子はいらない
        "movement": False,
        "velocity": "floating",
        "mass": 156.96/9.81,
        "COM": [0, 0, 0.25],
        "MOI": [50, 50, 50]
    }
    print('------------------------------------')
    for key, value in data.items():
        print(key, '\t', value)
    print('------------------------------------')
    f = open("./floatingbody.json", 'w')
    json.dump(data, f, ensure_ascii=True, indent=4)
    f.close()
except:
    print("no floatingbody")

# BEMの格子をCADで作る場合は，x,y,zそのままで作っていい．x,y,zがx,z,yで出力されたりはしない．
