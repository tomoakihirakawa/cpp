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
# output_directory = home+"/BEM/BEM_dt0d01_wavemake_Xdir_H0d10_L0d25_A0d01_ww1_1d0_using_Dombre2019_respect_Dirichlet_on_CORNER_EMT_K0d5_without_remeshing"

settingBEM = {
    "beta": 0.,
    "K": 1.,
    "max_dt": 0.02,
    "mesh": 24,
    "stop_remesh_time": 100.,
    "force_remesh_time": 0.1
}


mesh = settingBEM["mesh"]

output_directory = home + "/BEM/Li2002_h0_0d3048_shorttank" \
                        + "_K" + str(settingBEM["K"])\
                        + "_beta" + str(settingBEM["beta"])\
                        + "_max_dt" + str(settingBEM["max_dt"])\
                        + "_Mesh" + str(settingBEM["mesh"]) \
                        + "_stopremesh" + str(settingBEM["stop_remesh_time"])


folder = "../../obj/tank/Li2002_h0_0d3048_shorttank"


data = {
    "name": "water",
    "objfile": folder+"/water"+str(mesh)+".obj",
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

data = {
    "name": "wavemaker",
    "objfile": folder+"/wavemaker.obj",
    "radius": 0.01,
    "reverseNormal": True,
    "ignore": False,
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
    "objfile": folder+"/wavetank.obj",
    "radius": 0.01,
    "reverseNormal": True,
    "ignore": False,
    # "scale": [0.001, 0.001, 0.001, 0.001],  # モデルがmm単位なのでメートルに変換,
    # -------------------------------------------------------- #
    "output_vtu_file_name": "wavetank",  # 拡張子はいらない
    "output_pvd_file_name": "wavetank",  # 拡張子はいらない
    "output_directory": output_directory,
}

f = open("./wavetank.json", 'w')
json.dump(data, f, ensure_ascii=True, indent=4)
f.close()
