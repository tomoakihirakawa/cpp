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

output_directory = home+"/BEM/BEM_amp035_Dombre2019_move1"

data = {
    "name": "water",
    "objfile": "../../obj/tank/sawai/water_sawai_250x250mm15.obj",
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

data = {
    "name": "tank",
    "objfile": "../../obj/tank/sawai/tank_sawai_254x254mm15.obj",
    "radius": 0.01,
    # "scale": [0.001, 0.001, 0.001, 0.001],  # モデルがmm単位なのでメートルに変換,
    # -------------------------------------------------------- #
    "output_vtu_file_name": "tank",  # 拡張子はいらない
    "output_pvd_file_name": "tank",  # 拡張子はいらない
    "output_directory": output_directory,
}

f = open("./tank.json", 'w')
json.dump(data, f, ensure_ascii=True, indent=4)
f.close()
