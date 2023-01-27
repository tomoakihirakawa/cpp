
from math import pi
import json
import math
from os.path import expanduser
home = expanduser("~")
rho = 1000.
g = 9.81

'''
プログラムを回す際に面倒な事は，入力ファイルの設定方法．
入力ファイルの作り方をドキュメントで示されても，具体的な例がないとわかりにくい．
例があっても，例と違う場合どうすればいいかなど，わからないことは多い．
このように，入力ファイルを生成するプログラムを作っておけば，その面倒をだいぶ解消できる．
'''

#! -------------------------------------------------------- #

water = {
    "name": "water",
    "type": "Fluid",
    "output_vtu_file_name": "water",  # 拡張子はいらない
    "output_pvd_file_name": "water",  # 拡張子はいらない
}
wavemaker = {
    "name": "wavemaker",
    "type": "RigidBody",
    "output_vtu_file_name": "wavemaker",  # 拡張子はいらない
    "output_pvd_file_name": "wavemaker",  # 拡張子はいらない
    # "velocity": ["Goring1979", 4.],
    # "velocity": ["Retzler2000", .1],
    "velocity": ["Chaplin2000", .1],
}
wavetank = {
    "name": "wavetank",
    "type": "RigidBody",
    # "ignore": wavetank_ignore,
    "output_vtu_file_name": "wavetank",  # 拡張子はいらない
    "output_pvd_file_name": "wavetank",  # 拡張子はいらない
}
floatingbody = {
    "name": "floatingbody",
    "type": "RigidBody",
    "output_vtu_file_name": "floatingbody",  # 拡張子はいらない
    "output_pvd_file_name": "floatingbody",  # 拡張子はいらない
    "velocity": "floating",
}
floatingbodies = []
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

SimulationCase = "Retzler2000simple"
match SimulationCase:
    case "tsukada_wave_generation_test":
        settingBEM["inputfiles"] = ["water.json",
                                    "wavemaker.json", "wavetank.json"]
        settingBEM["max_dt"] = 0.2
        settingBEM["output_directory"] = home + "/BEM/2022Tsukada_A1_T5_300"
        #
        objfolder = "../../obj/2022Tsukada/wave_generation_test"
        water["objfile"] = objfolder+"/"+"water200.obj"
        wavemaker["objfile"] = objfolder+"/"+"paddle100.obj"
        wavetank["objfile"] = objfolder+"/"+"tank100.obj"
        wavetank_ignore = False
        start = 0
        A = 7.
        T = 8
        h = 110
        l = 0
        wavemaker["velocity"] = ["flap", start, A, T, h, l, 0, 1, 0]
    case "wave_generation_test":
        settingBEM["inputfiles"] = ["water.json",
                                    "wavemaker.json", "wavetank.json"]
        settingBEM["max_dt"] = 0.2
        settingBEM["output_directory"] = home + "/BEM/wave_generation_test"
        #
        # objfolder = "wave_generation_test"
        objfolder = "2022Tonegawa/test20221020"
        water["objfile"] = "../../obj/"+objfolder+"/"+"water200.obj"
        wavemaker["objfile"] = "../../obj/"+objfolder+"/"+"wavemaker0.obj"
        wavetank["objfile"] = "../../obj/"+objfolder+"/"+"tank0.obj"
        wavetank_ignore = False
        start = 0
        A = 7
        T = 8
        h = 150
        l = 0
        wavemaker["velocity"] = ["flap", start, A, T, h, l, 0, 1, 0]
    case "FloatingWind1":
        settingBEM["max_dt"] = 0.05
        settingBEM["output_directory"] = home + "/BEM/FloatingWind1"
        settingBEM["inputfiles"] = ["water.json",
                                    "floatingbody.json",
                                    "wavemaker.json"]
        #
        objfolder = "FloatingWind1"
        water["objfile"] = "../../obj/"+objfolder+"/"+"/water150_remeshed.obj"
        floatingbody["objfile"] = "../../obj/" + \
            objfolder+"/"+"floating_body.obj"
        floatingbody["mass"] = (0.5*0.5-0.3*0.3)*rho*0.1
        floatingbody["COM"] = [0., 0, 0.25]
        floatingbody["MOI"] = [20, 20, 20]
        wavemaker["objfile"] = "../../obj/"+objfolder+"/"+"tank.obj"
        wavemaker["velocity"] = ["Chaplin2000", 1.]
    case "FloatingWind2":
        settingBEM["max_dt"] = 0.05
        settingBEM["output_directory"] = home + "/BEM/FloatingWind2"
        settingBEM["inputfiles"] = ["water.json",
                                    "floatingbody.json",
                                    "wavemaker.json"]
        objfolder = "FloatingWind2"
        water["objfile"] = "../../obj/"+objfolder+"/"+"/water100.obj"
        floatingbody["objfile"] = "../../obj/" + \
            objfolder+"/"+"floatingbody60.obj"
        floatingbody["mass"] = (0.4*0.4 - 0.2*0.2)*rho*0.05
        floatingbody["COM"] = [0., 0., 0.]
        floatingbody["MOI"] = [20, 20, 20]
        wavemaker["objfile"] = "../../obj/"+objfolder+"/"+"tank.obj"
        wavemaker["velocity"] = ["Chaplin2000", 10.]
    case "Tsukada2022":
        settingBEM["max_dt"] = 0.5
        settingBEM["output_directory"] = home + "/BEM/2022Tsukada"
        settingBEM["inputfiles"] = ["water.json",
                                    "floatingbody.json",
                                    "wavemaker.json"]
        objfolder = "2022Tsukada"
        water["objfile"] = "../../obj/"+objfolder+"/"+"/water270.obj"
        floatingbody["objfile"] = "../../obj/" + \
            objfolder+"/"+"floatingbody1.obj"
        floatingbody["mass"] = m = 8.56026e+07/g
        floatingbody["COM"] = [0., 0., 52.5+7.275]
        floatingbody["radius of gyration"] = [18.7, 20.3, 21.1]
        floatingbody["MOI"] = [m*math.pow(floatingbody["radius of gyration"][0], 2),
                               m*math.pow(floatingbody["radius of gyration"][0], 2),
                               m*math.pow(floatingbody["radius of gyration"][0], 2)]
        wavemaker["objfile"] = "../../obj/"+objfolder+"/"+"tank.obj"
        wavemaker["velocity"] = ["Chaplin2000", 0.1, .5, 2*math.pi/5.]
    case "Tonegawa2022":
        settingBEM["max_dt"] = 0.05
        settingBEM["output_directory"] = home + "/BEM/2022Tonegawa"
        settingBEM["inputfiles"] = ["water.json",
                                    "floatingbody0.json",
                                    "floatingbody1.json",
                                    "floatingbody2.json",
                                    "floatingbody3.json",
                                    "wavetank.json",
                                    "wavemaker.json"]
        objfolder = "2022Tonegawa"
        water["objfile"] = "../../obj/"+objfolder+"/water90.obj"
        for i in range(4):
            m = (0.5*0.5)*0.1*rho
            radius_of_gyration = [100., 100., 100.]
            name = "floatingbody" + str(i)
            floatingbodies.append(
                {
                    "name": name,
                    "type": "RigidBody",
                    "output_vtu_file_name": name,  # 拡張子はいらない
                    "output_pvd_file_name": name,  # 拡張子はいらない
                    "velocity": "fixed",
                    "objfile": "../../obj/"+objfolder+"/"+name+".obj",
                    "mass": m,
                    "expected_force": m*g,
                    "COM": [-3+1.5*(i % 5), 2.6-2.6*int(i / 5),  2],
                    "MOI": [m*math.pow(radius_of_gyration[0], 2),
                            m*math.pow(radius_of_gyration[1], 2),
                            m*math.pow(radius_of_gyration[2], 2)]
                })
        wavemaker["objfile"] = "../../obj/"+objfolder+"/"+"wavemaker.obj"
        wavemaker["velocity"] = ["Chaplin2000", 0.1, 0.05, 2*math.pi/1.]
        wavetank["objfile"] = "../../obj/"+objfolder+"/"+"wavetank.obj"
    case "MultipleFloatingBodies":
        # b!=========================================================================
        settingBEM["max_dt"] = 0.05
        settingBEM["output_directory"] = home + "/BEM/MultipleFloatingBodies"
        settingBEM["inputfiles"] = ["water.json",
                                    "wavemaker.json",
                                    "wavetank.json",
                                    "floatingbody1.json",
                                    "floatingbody2.json",
                                    "floatingbody3.json",
                                    "floatingbody4.json",
                                    "floatingbody5.json",
                                    "floatingbody6.json",
                                    "floatingbody7.json",
                                    "floatingbody8.json",
                                    "floatingbody9.json"]
        # -------------------------------------------------------- #
        objfolder = "MultipleFloatingBodies"
        water["objfile"] = "../../obj/"+objfolder+"/"+"/water300remeshed.obj"
        k = 0
        for j in range(1, 4):
            for i in range(1, 4):
                k = 1+k
                m = (0.5*0.5)*0.3*rho
                radius_of_gyration = [100., 100., 100.]
                name = "floatingbody" + str(k)
                floatingbodies.append(
                    {
                        "name": name,
                        "type": "RigidBody",
                        "output_vtu_file_name": name+"_",  # 拡張子はいらない
                        "output_pvd_file_name": name+"_",  # 拡張子はいらない
                        "velocity": "floating",
                        "objfile": "../../obj/"+objfolder+"/"+name+".obj",
                        "mass": m,
                        "expected_force": m*g,
                        "COM": [-5+5*(i-1), 4-4*(j-1),  2],
                        "MOI": [m*math.pow(radius_of_gyration[0], 2),
                                m*math.pow(radius_of_gyration[1], 2),
                                m*math.pow(radius_of_gyration[2], 2)]
                    })
        wavemaker["objfile"] = "../../obj/"+objfolder+"/"+"wavemaker.obj"
        wavemaker["velocity"] = ["Chaplin2000", 0.1, 0.05, 2*math.pi/1.]
        wavetank["objfile"] = "../../obj/"+objfolder+"/"+"tank.obj"
    case "sloshing_tank_torus":
        # water = "/water"+str(settingBEM["mesh"])+".obj"
        water = "/water100.obj"
        objfolder = "sloshing_tank_torus"
        wavemaker = "torus4.obj"
        wavetank = "tank3.obj"
        wavetank_ignore = False
        settingBEM["output_directory"] = home + "/BEM/tourus"
        settingBEM["inputfiles"] = ["water.json",
                                    "floatingbody.json", "wavetank.json"]
    case "sloshing_tank_sphere":
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
    case "sloshing_tank":
        settingBEM["max_dt"] = 0.01
        objfolder = "sloshing_tank"
        settingBEM["output_directory"] = home + "/BEM/sloshing_tank"
        settingBEM["inputfiles"] = ["water.json",
                                    "wavetank.json"]
        water["objfile"] = "../../obj/"+objfolder+"/"+"/water20.obj"
        wavetank["objfile"] = "../../obj/"+objfolder+"/"+"tank3.obj"
        wavetank["velocity"] = ["Chaplin2000", 0.1, 0.05, 2*math.pi]
    case "Retzler2000simple":
        # objfolder = "chaplin2000simple"
        # water = "/water.obj"
        # wavemaker = "cylinder100.obj"
        # wavetank = "tank.obj"
        # wavetank_ignore = False
        # settingBEM["inputfiles"] = ["water.json",
        #                             "wavemaker.json", "wavetank.json"]
        # settingBEM["max_dt"] = 0.002
        # settingBEM["output_directory"] = home + "/BEM/Retzler2000_" + \
        #     str(settingBEM["max_dt"]).replace(
        #         '.', 'd') + "AreaWeightedVelocity_1_Linear_simple"
        # output_directory = home + "/BEM/Retzler2000_0d005_linear"
        # output_directory = home + "/BEM/test"

        settingBEM["inputfiles"] = ["water.json",
                                    "wavemaker.json", "wavetank.json"]
        settingBEM["max_dt"] = 0.002
        settingBEM["output_directory"] = home + "/BEM/Retzler2000_" + \
            str(settingBEM["max_dt"]).replace('.', 'd') + \
            "_Linear_10times_20221016"
        #
        objfolder = "chaplin2000simple"
        water["objfile"] = "../../obj/"+objfolder+"/"+"water.obj"
        wavemaker["objfile"] = "../../obj/"+objfolder+"/"+"cylinder100.obj"
        wavetank["objfile"] = "../../obj/"+objfolder+"/"+"tank.obj"
        wavetank_ignore = False
        wavemaker["velocity"] = ["Retzler2000", 0.1]
        # output_directory = home + "/BEM/Retzler2000_0d005_linear"
        # output_directory = home + "/BEM/test"
    case "Goring1979":
        settingBEM["max_dt"] = 0.04
        settingBEM["inputfiles"] = ["water.json",
                                    "wavemaker.json", "wavetank.json"]
        mesh = "250"
        settingBEM["output_directory"] = home + "/BEM/Goring1979_" + \
            str(settingBEM["max_dt"]).replace('.', 'd') + \
            "_Linear" + mesh
        #
        objfolder = "Goring1979"
        water["objfile"] = "../../obj/" + \
            objfolder+"/"+"/water_short"+mesh+".obj"
        wavemaker["objfile"] = "../../obj/"+objfolder+"/"+"wavemaker_short.obj"
        wavemaker["velocity"] = ["Goring1979", 4.]
        wavetank["objfile"] = "../../obj/"+objfolder+"/"+"wavetank_short.obj"
    case "Li2002":
        water = "/water60_refined.obj"
        objfolder = "/Li2002"
        wavemaker = "wavemaker.obj"
        wavetank = "tank.obj"
        wavetank_ignore = False
        settingBEM["inputfiles"] = ["water.json",
                                    "wavemaker.json", "wavetank.json"]
        settingBEM["max_dt"] = 0.01
        settingBEM["output_directory"] = home + "/BEM/Li2002_" + str(settingBEM["max_dt"]).replace('.', 'd') + \
            "_Hh0d2_AreaWeightedVelocity_1_Linear1"
        # output_directory = home + "/BEM/Retzler2000_0d005_linear"
        # output_directory = home + "/BEM/test"
    case "Retzler2000":
        settingBEM["inputfiles"] = ["water.json",
                                    "wavemaker.json", "wavetank.json"]
        settingBEM["max_dt"] = 0.002
        settingBEM["output_directory"] = home + "/BEM/Retzler2000_" + \
            str(settingBEM["max_dt"]).replace('.', 'd') + \
            "_Linear_10times_"
        #
        objfolder = "chaplin2000"
        water["objfile"] = "../../obj/"+objfolder+"/"+"water_remeshed2.obj"
        wavemaker["objfile"] = "../../obj/"+objfolder+"/"+"cylinder100.obj"
        wavetank["objfile"] = "../../obj/"+objfolder+"/"+"tank3.obj"
        wavetank_ignore = False
        wavemaker["velocity"] = ["Retzler2000", 0.1]
        # output_directory = home + "/BEM/Retzler2000_0d005_linear"
        # output_directory = home + "/BEM/test"
    case "Chaplin2000":
        settingBEM["inputfiles"] = ["water.json",
                                    "wavemaker.json", "wavetank.json"]
        settingBEM["max_dt"] = 0.01
        settingBEM["output_directory"] = home + "/BEM/Chaplin2000_" + \
            str(settingBEM["max_dt"]).replace('.', 'd') + "_Linear_10times_"
        objfolder = "chaplin2000"
        water["objfile"] = "../../obj/"+objfolder+"/"+"water_remeshed2.obj"
        wavemaker["objfile"] = "../../obj/"+objfolder+"/"+"cylinder100.obj"
        wavetank["objfile"] = "../../obj/"+objfolder+"/"+"tank3.obj"
        wavetank_ignore = False

        h = 0.5
        w = 1.257 * math.sqrt(g/h)
        # Hu (2002) shows w = 5.57065
        # A = 0.046 * h
        A = 0.08 * h
        print("w=", w, ", A=", A, ", T=", 2*math.pi/w)
        wavemaker["velocity"] = ["Chaplin2000", 0., A, w]
        # output_directory = home + "/BEM/Retzler2000_0d005_linear"
        # output_directory = home + "/BEM/test"
        # # water = "/water30.obj"
        # # water = "/water115_test.obj"
        # water = "/water_remeshed2.obj"
        # # water = "/water80_refined_120.obj"
        # # water = "/water80_test.obj"
        # # water = "/water100_test.obj"
        # # water = "/water100_refined_80.obj"
        # # water = "/water30_refined_199.obj"
        # # water = "/water30_refined_91.obj"
        # objfolder = "chaplin2000"
        # wavemaker = "cylinder100.obj"
        # wavetank = "tank3.obj"
        # wavetank_ignore = False
        # settingBEM["inputfiles"] = ["water.json",
        #                             "wavemaker.json", "wavetank.json"]
        # settingBEM["max_dt"] = 0.01
        # settingBEM["output_directory"] = home + "/BEM/Chaplin2000_" + \
        #     str(settingBEM["max_dt"]).replace('.', 'd') + \
        #     "_Linear_10times"
        # # output_directory = home + "/BEM/Retzler2000_0d005_linear"
        # # output_directory = home + "/BEM/test"
    case "Chaplin2000simple":
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
            str(settingBEM["max_dt"]).replace('.', 'd') + "_linear3_"
    case "refinement":
        water = "/water_remeshed2.obj"
        objfolder = "chaplin2000"
        wavemaker = "cylinder100.obj"
        wavetank = "tank3.obj"
        wavetank_ignore = False
        settingBEM["output_directory"] = home + "/BEM/refinement"
        settingBEM["grid_refinement"] = 200
        settingBEM["inputfiles"] = ["water.json",
                                    "wavemaker.json", "wavetank.json"]

#@ -------------------------------------------------------- #
#@                  setting.json を出力                      #
#@ -------------------------------------------------------- #
print('------------------------------------')
for key, value in settingBEM.items():
    print(f'{key: <{20}}', '\t', '\033[91m', value, '\033[0m')
print('------------------------------------')
f = open("./inputs/setting.json", 'w')
json.dump(settingBEM, f, ensure_ascii=True, indent=4)
f.close()

#@ -------------------------------------------------------- #
#@           その他，water.json,tank.json などを出力           #
#@ -------------------------------------------------------- #
for INPUTS in [water, wavemaker, wavetank, floatingbody]:
    if INPUTS["name"]+".json" in settingBEM["inputfiles"]:
        print('------------------------------------')
        for key, value in INPUTS.items():
            print(f'{key: <{20}}', '\t', '\033[92m', value, '\033[0m')
        print('------------------------------------')
        f = open("./inputs/"+INPUTS["name"]+".json", 'w')
        json.dump(INPUTS, f, ensure_ascii=True, indent=4)
        f.close()

for INPUTS in floatingbodies:
    if INPUTS["name"]+".json" in settingBEM["inputfiles"]:
        print('------------------------------------')
        for key, value in INPUTS.items():
            print(f'{key: <{20}}', '\t', '\033[92m', value, '\033[0m')
        print('------------------------------------')
        f = open("./inputs/"+INPUTS["name"]+".json", 'w')
        json.dump(INPUTS, f, ensure_ascii=True, indent=4)
        f.close()

# BEMの格子をCADで作る場合は，x,y,zそのままで作っていい．x,y,zがx,z,yで出力されたりはしない．
