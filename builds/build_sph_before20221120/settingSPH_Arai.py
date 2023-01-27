from math import pi
import json
import math

# edge_length = 2.*0.05

density = 1000.
graity = 9.81

H = 200/1000
# particle_spacing = 0.008
particle_spacing = .5
data = {
    "density": density,
    # -------------------------------------------------------- #
    "initial_surface_z_position": H,
    # ---------------- 時間間隔dtに関する設定値 ----------------- #
    # 最初のstep==0の場合はdt=1E-10とする
    "C_CFL_velocity": 0.075,  # dt = C_CFL_velocity*h/Max(U)
    "C_CFL_accel": 0.25,     # dt = C_CFL_accel*sqrt(h/Max(A))
    # ------------------ 粒子配置に関する設定値 ------------------ #
    "particle_spacing": particle_spacing,
    # -------------------------------------------------------- #
    # バケットのバウンディングボックスの外に達した流体粒子は削除する．
    "buckets_xbounds": [-50., 50.],
    "buckets_ybounds": [-50., 50.],
    "buckets_zbounds": [-10., 10.],
    #@ ---------------------- 平滑化半径に関するの設定値（計算精度に関わる） --------------------- #
    "C_SML_sigma":3.,  # 一般的な平滑化距離．5次のスプラインの場合3h離れた粒子は影しない
    "C_SML_h": 1.,  # 一般的な平滑化距離．5次のスプラインの場合3h離れた粒子は影しない
    # "C_SML": 3.*.9,  # 一般的な平滑化距離．5次のスプラインの場合3h離れた粒子は影しない
    "kNS_SML": 5.,  # k-nearest search. dxを決めるための近傍粒子数
    # radius_SPHはC_SML*dx
    # -------------------------------------------------------- #
    "mu": 0.001005,
    # -------------------------------------------------------- #
    "C_Tait": 10.*math.sqrt(2.*graity*H),  # テイトの式
    # ---------------------- 人口粘性係数 ---------------------- 
    "C_artificial_viscousity_alpha": 0.001,
    "C_artificial_viscousity_beta": 0.,
    # ------------------------ 準備時間 ------------------------ #
    "preparation_max_dt": 0.0005,
    "preparation_time": 1.,
    "preparation_time_step": 200,
    "preparation_C_artificial_viscousity_alpha": 0.03,
    "preparation_C_artificial_viscousity_beta": 0,
    # -------------------------------------------------------- #
    "inputfiles" : ["water.json", "tank.json"]
}

data["max_dt"] = data["C_CFL_velocity"]*data["C_SML_h"]*particle_spacing/(data["C_Tait"]/10.)
data["max_dt"] = 0.0005

print('------------------------------------')
for key, value in data.items():
    print(key,'\t',value)
print('------------------------------------')
f = open("./settingSPH.json", 'w')
json.dump(data, f, ensure_ascii=True, indent=4)
f.close()

#! -------------------------------------------------------- #
#! -------------------------------------------------------- #
#! -------------------------------------------------------- #

data = {
    "name": "water",    
    "objfile": "../../obj/2022Arai/water_simple.obj",
}
print('------------------------------------')
for key, value in data.items():
    print(f'{key: <{20}}', '\t', '\033[91m', value, '\033[0m')
print('------------------------------------')

f = open("./water.json", 'w')
json.dump(data, f, ensure_ascii=True, indent=4)
f.close()

#! -------------------------------------------------------- #
#! -------------------------------------------------------- #
#! -------------------------------------------------------- #

data = {
    "name": "wave_maker",    
    "objfile": "../../obj/2022Arai/tank70.obj",
    # "rotate": [0, 1, 0, 0],
    # "reverseNormal": True,
    # "scale": [0, 0, 0, 0],  # モデルがmm単位なのでメートルに変換,
    # "type": "RigidBody",
    "depth_list": [-particle_spacing/2.,
                   -particle_spacing/2.*3.],
    # "volume_of_a_particle": volume_of_a_particle,
    "ignore":False,
    "density": density
}
print('------------------------------------')
for key, value in data.items():
    print(f'{key: <{20}}', '\t', '\033[92m', value, '\033[0m')
print('------------------------------------')

f = open("./wave_maker.json", 'w')
json.dump(data, f, ensure_ascii=True, indent=4)
f.close()

#! -------------------------------------------------------- #
#! -------------------------------------------------------- #
#! -------------------------------------------------------- #

data = {
    "name": "tank",
    "objfile": "../../obj/2022Arai/tank10.obj",
    "type": "RigidBody",
    "depth_list": [-particle_spacing/2.,
                   -particle_spacing/2.*3.],
    # "volume_of_a_particle": volume_of_a_particle,
    "ignore":False,
    "density": density
}
print('------------------------------------')
for key, value in data.items():
    print(f'{key: <{20}}', '\t', '\033[92m', value, '\033[0m')
print('------------------------------------')

f = open("./tank.json", 'w')
json.dump(data, f, ensure_ascii=True, indent=4)
f.close()


#! -------------------------------------------------------- #
#! -------------------------------------------------------- #
#! -------------------------------------------------------- #

data = {
    "name": "tank",
    # "objfile": "../../obj/cube0.obj",
    "objfile": "../../obj/2022Arai/tank.obj",
    # "rotate": [math.pi/3., 1, 1, 0],
    # "translate": [1., 1., 1.],
    # "reverseNormal": True,
    "scale": [1./1.5, 0,0,0],  # モデルがmm単位なのでメートルに変換,
    "translate": [0, 0, -0.3],
    "ignore":True,
    "type": "RigidBody",
    "depth_list": [-particle_spacing/2.],
    # "volume_of_a_particle": volume_of_a_particle,
    "density": density

}
print('------------------------------------')
for key, value in data.items():
    print(f'{key: <{20}}', '\t', '\033[92m', value, '\033[0m')
print('------------------------------------')
f = open("./tank_init.json", 'w')
json.dump(data, f, ensure_ascii=True, indent=4)
f.close()


