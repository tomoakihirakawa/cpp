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
output_directory = home + \
    "/BEM/sloshing_with_pole/BEM_dt0d01_Xdir_H0d10_L0d25_A0d01_ww1_1d0_EMT_K0d5"
# output_directory = home+"/BEM/test"
data = {
    "name": "water",
    "objfile": "../../obj/tank/sloshing_with_pole/water14.obj",
    "radius": 0.01,
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
    "objfile": "../../obj/tank/sloshing_with_pole/tank13.obj",
    "ignore": False,
    "radius": 0.01,
    "center_of_mass": [0, 0, 0],
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
    "objfile": "../../obj/tank/sloshing_with_pole/tank.obj",
    "ignore": True,
    "radius": 0.01,
    "reverseNormal": True,
    "center_of_mass": [0, 0, 0],
    # "scale": [0.001, 0.001, 0.001, 0.001],  # モデルがmm単位なのでメートルに変換,
    # -------------------------------------------------------- #
    "output_vtu_file_name": "wavetank",  # 拡張子はいらない
    "output_pvd_file_name": "wavetank",  # 拡張子はいらない
    "output_directory": output_directory,
}

f = open("./wavetank.json", 'w')
json.dump(data, f, ensure_ascii=True, indent=4)
f.close()
