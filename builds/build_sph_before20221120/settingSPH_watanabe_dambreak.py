from math import pi
import json
import math

# edge_length = 2.*0.05

density = 1000.
graity = 9.81

H = 292/1000
particle_spacing = 0.008
data = {
    "density": density,
    # -------------------------------------------------------- #
    "initial_surface_z_position": H + 500/1000,
    # ---------------- 時間間隔dtに関する設定値 ----------------- #
    # 最初のstep==0の場合はdt=1E-10とする
    "C_CFL_velocity": 0.075,  # dt = C_CFL_velocity*h/Max(U)
    "C_CFL_accel": 0.25,     # dt = C_CFL_accel*sqrt(h/Max(A))
    # ------------------ 粒子配置に関する設定値 ------------------ #
    "particle_spacing": particle_spacing,
    "xbounds": [659/1000 + particle_spacing/2, 805/1000 - particle_spacing/2],
    "ybounds": [-75/1000 + particle_spacing/2, 75/1000 - particle_spacing/2],
    "zbounds": [500/1000 + particle_spacing/2, 500/1000+H + particle_spacing/2],
    # -------------------------------------------------------- #
    # バケットのバウンディングボックスの外に達した流体粒子は削除する．
    "buckets_xbounds": [0., 0.9],
    "buckets_ybounds": [-.14, .14],
    "buckets_zbounds": [0.4, 1.2],
    #@ ---------------------- 平滑化半径に関するの設定値（計算精度に関わる） --------------------- #
    "C_SML_sigma":2.,  # 一般的な平滑化距離．5次のスプラインの場合3h離れた粒子は影しない
    "C_SML_h": 1.09,  # 一般的な平滑化距離．5次のスプラインの場合3h離れた粒子は影しない
    # "C_SML": 3.*.9,  # 一般的な平滑化距離．5次のスプラインの場合3h離れた粒子は影しない
    "kNS_SML": 5.,  # k-nearest search. dxを決めるための近傍粒子数
    # -------------------------------------------------------- #
    "mu": 0.001005,
    # -------------------------------------------------------- #
    "C_Tait": 10.*math.sqrt(2.*9.81*H),  # テイトの式
    # ---------------------- 人口粘性係数 ---------------------- 
    "C_artificial_viscousity_alpha": 0.001,
    "C_artificial_viscousity_beta": 0.,
    # ------------------------ 準備時間 ------------------------ #
    "preparation_max_dt": 0.00005,
    "preparation_time": 1.,
    "preparation_time_step": 200,
    "preparation_C_artificial_viscousity_alpha": 0.03,
    "preparation_C_artificial_viscousity_beta": 0,
    # -------------------------------------------------------- #
}

data["max_dt"] = data["C_CFL_velocity"]*data["C_SML_h"]*particle_spacing/(data["C_Tait"]/10.)

total_volume = 1
total_volume *= abs(data["xbounds"][1]-data["xbounds"][0]+ particle_spacing)
total_volume *= abs(data["ybounds"][1]-data["ybounds"][0]+ particle_spacing)
total_volume *= abs(data["zbounds"][1]-data["zbounds"][0]+ particle_spacing)
print("体積", total_volume)
data.update({"total_volume": total_volume})

xlen = data["xbounds"][1]-data["xbounds"][0]
ylen = data["ybounds"][1]-data["ybounds"][0]
zlen = data["zbounds"][1]-data["zbounds"][0]
dx = data["particle_spacing"]
total_particle = (round(xlen/dx+1)*round(ylen/dx+1)*round(zlen/dx+1))
volume_of_a_particle = total_volume/total_particle
data.update({"volume_of_a_particle": total_particle})
print("粒子点数", total_particle)
print("各粒子体積", volume_of_a_particle)
print('------------------------------------')
for key, value in data.items():
    print(key,'\t',value)
print('------------------------------------')
f = open("./settingSPH.json", 'w')
json.dump(data, f, ensure_ascii=True, indent=4)
f.close()

print("rho*g*h = ", 1000.*9.81*H)

#! -------------------------------------------------------- #
#! -------------------------------------------------------- #
#! -------------------------------------------------------- #

data = {
    "name": "wave_maker",
    "objfile": "../../obj/tank/sloshing_tank.obj",
    # "rotate": [0, 1, 0, 0],
    "reverseNormal": True,

    "scale": [1/100, 0, 0, 0],  # モデルがmm単位なのでメートルに変換,
    "translate": [0, 0, 2.5],  # モデルがmm単位なのでメートルに変換,

    "depth_list": [-particle_spacing/2.],
    "volume_of_a_particle": volume_of_a_particle,
    "density": density
}

f = open("./wave_maker.json", 'w')
json.dump(data, f, ensure_ascii=True, indent=4)
f.close()

#! -------------------------------------------------------- #
#! -------------------------------------------------------- #
#! -------------------------------------------------------- #

data = {
    "name": "tank",
    "objfile": "../../obj/watanabe2021/tank_for_dambreak_experiment_Koshizuka_Moyce15.obj",
    # "rotate": [2*math.pi, 1, 0, 0],
    "scale": [1/1000, 0, 0, 0],  # モデルがmm単位なのでメートルに変換,
    "depth_list": [-particle_spacing/2.,
                   -particle_spacing/2.*3],
    "volume_of_a_particle": volume_of_a_particle,
    "density": density
}

f = open("./tank.json", 'w')
json.dump(data, f, ensure_ascii=True, indent=4)
f.close()


#! -------------------------------------------------------- #
#! -------------------------------------------------------- #
#! -------------------------------------------------------- #

data = {
    "name": "tank",
    "objfile": "../../obj/watanabe2021/tank_for_dambreak_experiment_Koshizuka_Moyce_init15.obj",
    # "rotate": [math.pi/2., 1, 0, 0],
    "scale": [1/1000, 0, 0, 0],  # モデルがmm単位なのでメートルに変換,
    "depth_list": [-particle_spacing/2.,
                   -particle_spacing/2.*3],
    "volume_of_a_particle": volume_of_a_particle,
    "density": density

}

f = open("./tank_init.json", 'w')
json.dump(data, f, ensure_ascii=True, indent=4)
f.close()
