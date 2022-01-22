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
    "beta": 0.,  # 拡張子はいらない
    "K": 1.,  # 拡張子はいらない
    "max_dt": 0.02,  # 拡張子はいらない
    "mesh": 24,
    "stop_remesh_time": 100.
}

output_directory = home + "/BEM/XueAndLin2011_large_amplitude" \
                        + "_K" + str(settingBEM["K"])\
                        + "_beta" + str(settingBEM["beta"])\
                        + "_max_dt" + str(settingBEM["max_dt"])\
                        + "_Mesh" + str(settingBEM["mesh"]) \
                        + "_stopremesh" + str(settingBEM["stop_remesh_time"]) \
                        + "_not_respect_Dirichlet_on_CORNER"

settingBEM["output_directory"] = output_directory


f = open("./settingBEM.json", 'w')
json.dump(settingBEM, f, ensure_ascii=True, indent=4)
f.close()


dir = "../../obj/XueAndLin2011/"
#! -------------------------------------------------------- #
data = {
    "name": "water",
    "objfile": dir + "water"+str(settingBEM["mesh"])+".obj",
    "radius": 0.01,
    # "rotate": [math.pi/2., 1, 0, 0],
    # "scale": [0.001, 0.001, 0.001, 0.001],  # モデルがmm単位なのでメートルに変換,
    # -------------------------------------------------------- #
    "output_vtu_file_name": "water",  # 拡張子はいらない
    "output_pvd_file_name": "water",  # 拡張子はいらない
    "output_directory": output_directory,
}

f = open("./water.json", 'w')
json.dump(data, f, ensure_ascii=True, indent=4)
f.close()

#! -------------------------------------------------------- #

# data = {
#     "name": "tank",
#     "objfile": "../../obj/tank/sawai/tank_sawai_254x254mm5.obj",
#     "radius": 0.01,
#     # "scale": [0.001, 0.001, 0.001, 0.001],  # モデルがmm単位なのでメートルに変換,
#     # -------------------------------------------------------- #
#     "output_vtu_file_name": "tank",  # 拡張子はいらない
#     "output_pvd_file_name": "tank",  # 拡張子はいらない
#     "output_directory": output_directory,
# }

# f = open("./tank.json", 'w')
# json.dump(data, f, ensure_ascii=True, indent=4)
# f.close()


#! -------------------------------------------------------- #

data = {
    "name": "wavemaker",
    "objfile": dir + "tank1.obj",
    "radius": 0.01,
    # "rotate": [math.pi/2., 1, 0, 0],
    "ignore": False,
    # "reverseNormal": True,
    # "scale": [0.001, 0.001, 0.001, 0.001],  # モデルがmm単位なのでメートルに変換,
    # -------------------------------------------------------- #
    "output_vtu_file_name": "wavemaker",  # 拡張子はいらない
    "output_pvd_file_name": "wavemaker",  # 拡張子はいらない
    "output_directory": output_directory,
}

f = open("./wavemaker.json", 'w')
json.dump(data, f, ensure_ascii=True, indent=4)
f.close()

#! -------------------------------------------------------- #

data = {
    "name": "wavetank",
    "objfile": dir + "tank.obj",
    "radius": 0.01,
    "ignore": True,
    "reverseNormal": True,
    # "rotate": [math.pi/2., 1, 0, 0],
    # "scale": [0.001, 0.001, 0.001, 0.001],  # モデルがmm単位なのでメートルに変換,
    # -------------------------------------------------------- #
    "output_vtu_file_name": "wavetank",  # 拡張子はいらない
    "output_pvd_file_name": "wavetank",  # 拡張子はいらない
    "output_directory": output_directory,
}

f = open("./wavetank.json", 'w')
json.dump(data, f, ensure_ascii=True, indent=4)
f.close()
