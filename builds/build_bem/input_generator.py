r'''DOC_EXTRACT 2_0_0_input_generator

# Input Generator

## For Ubuntu

はじめに，以下のコマンドを実行して，必要なパッケージをインストールする．

```shell
python3.11 -m pip install numpy==1.25.0
```

## For Mac OS X

`input_generator.py`は，BEM-MELの入力ファイルを生成するためのスクリプトである．
`input_generator.py`を実行する際に，オプションを加えることで，シミュレーションケースやメッシュの名前，波を作る方法，要素の種類，時間刻み幅，シフトさせる面の補間方法，接尾語，波の高さ，出力ディレクトリ
などを指定することができるようにしている．
これは，`input_generator.py`を直接編集することなく，コマンドライン引数を使って計算条件を変更することができるようにするためで，
また，これによって，複数の計算条件を一括で実行することができるようになる．

```shell
python3.11 input_generator.py -h
```

例えば，次のようにシェルスクリプトを作成し，`input_generator.py`を実行することで，入力ファイルを生成することができる．

```shell
python3.11 input_generator.py Tanizawa1996 -m water_no_float0d08 -e linear -wavemaker flap -dt 0.03 -H 0.05 -ALE linear -ALEPERIOD 1 -outputdir ~/BEM
```

```shell
#!/bin/sh
outputdir=${HOME}/BEM
case=Palm2016
mesh_array=("water_mod")
dt_array=(0.05 0.03)
H_array=(0.05 0.1)

for dt in ${dt_array[@]};do
    for H in ${H_array[@]};do
        for mesh in ${mesh_array[@]};do
            python3.11 input_generator.py ${case} -m ${mesh} -e pseudo_quad -wavemaker flap -dt ${dt} -H ${H} -ALE linear -ALEPERIOD 1 -outputdir ${outputdir}
            python3.11 input_generator.py ${case} -m ${mesh} -e linear -wavemaker flap -dt ${dt} -H ${H} -ALE linear -ALEPERIOD 1 -outputdir ${outputdir}
        done
    done
done
```

python3.11 input_generator.py Tanizawa1996 -m water_no... -e pseudo_quad -wavemaker potential -dt 0.01 -ALE linear

'''

import copy
import platform
import json
import math
import os
import sys
import argparse

from math import pi
from os.path import expanduser

sys.path.append("..")

input_dir = "./input_files/"
sys_home_dir = expanduser("~")
current_dir = os.path.dirname(os.path.abspath(__file__))
code_home_dir = os.path.join(current_dir, '../../../')

# 
from IOgenerator import generate_input_files, read_obj, Norm3d, Differece3d, Add3d

rho = 1000.
g = 9.81

# ---------------------------------------------------------------------------- #

parser = argparse.ArgumentParser(description='コマンドライン引数を受け取るためのパーサー')
parser.add_argument('case', type=str, help='シミュレーションケース')
parser.add_argument('-m', '--mesh', type=str, help='メッシュの名前')
parser.add_argument('-wavemaker', type=str, help='波を作る方法')
parser.add_argument('-e', '--element', default='linear', type=str, help='要素の種類 (default: linear)')
parser.add_argument('-dt', '--max_dt', type=float, required=True, help='時間刻み幅')
parser.add_argument('-ALE', type=str, default='pseudo_quad', help='シフトさせる面の補間方法 (default: pseudo_quad)')
parser.add_argument('-ALEPERIOD', type=str, default='1', help='ALEの周期 (default: 1)')
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
T = None
padding = 20

SimulationCase = args.case
print(f"{'SimulationCase:':>{padding}} {SimulationCase}")
element = args.element
print(f"{'element:':>{padding}} {element}")
ALE = args.ALE
print(f"{'ALE:':>{padding}} {ALE}")
ALEPERIOD = args.ALEPERIOD
print(f"{'ALEPERIOD:':>{padding}} {ALEPERIOD}")

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


def IO_dir(id):
    input_dir = "./input_files/" + id
    os.makedirs(input_dir, exist_ok=True)
    if default_output_dir is None:
        output_dir = sys_home_dir + "/BEM/" + id
    else:
        output_dir = default_output_dir + "/" + id
    # output_dir = "/Volumes/home/BEM/benchmark20231206/" + id
    os.makedirs(output_dir, exist_ok=True)
    return input_dir, output_dir

# ---------------------------------------------------------------------------- #
# add id


def add_id(id):
    if dt is not None:
        max_dt = dt
    else:
        max_dt = 1/20
    id += "_DT" + str(max_dt).replace(".", "d")

    if meshname != "":
        id += "_MESH" + meshname

    if element is not None:
        id += "_ELEM" + element

    if ALE is not None:
        id += "_ALE" + ALE

    if ALEPERIOD != "":
        id += "_ALEPERIOD" + ALEPERIOD

    if suffix != "":
        id += "_" + suffix
    return id
# ---------------------------------------------------------------------------- #

if "looping" in SimulationCase:
    # objfolder = code_home_dir + "/cpp/obj/WaveGeneration"

    objfolder = code_home_dir + "/cpp/obj/looping"

    water = {"name": "water",
             "type": "Fluid",
             "objfile": objfolder + "/water_remesh.obj"}

    vertices, triangles = read_obj(water["objfile"])
    water["vertices"] = vertices
    water["triangles"] = triangles

    tank = {"name": "tank", "type":
            "RigidBody", "isFixed": True,
            "objfile": objfolder + "/tank.obj"}

    start = 0.

    T = 1.
    H = 0.05
    a = H/2  # Ren 2015 used H=[0.1(a=0.05), 0.03(a=0.06), 0.04(a=0.02)]
    h = 1.5

    # wavemaker_type = "piston"
    # wavemaker_type = "flap"
    # wavemaker_type = "potential"

    id = SimulationCase
    id += "_H" + str(H).replace(".", "d")
    id += "_T"+str(T).replace(".", "d")
    id += "_h"+str(h).replace(".", "d")

    z_surface = h

    if "piston" in SimulationCase:
        wavemaker = {"name": "wavemaker",
                     "type": "RigidBody",
                     "velocity": ["piston", start, a, T, h, 1, 0, 0],
                     "objfile": f"{objfolder}/wavemaker.obj"}

    elif "flap" in SimulationCase:
        wavemaker = {"name": "wavemaker",
                     "type": "RigidBody",
                     "velocity": ["flap", start, a, T, h, h, 0, 1, 0],
                     "objfile": f"{objfolder}/wavemaker.obj"}
    else:
        wavemaker = {"name": "wavemaker",
                     "type": "SoftBody",
                     "isFixed": True,
                     "velocity": ["velocity", start, a, T, h, z_surface],
                     "objfile": f"{objfolder}/wavemaker.obj"}

    inputfiles = [tank, wavemaker, water]

    setting = {"max_dt": 0.03,
               "end_time_step": 10000,
               "end_time": 9}

    generate_input_files(inputfiles, setting, IO_dir, id)
elif "Tanizawa1996" in SimulationCase:

    r'''DOC_EXTRACT 2_1_0_validation_Tanizawa1996

    ## Tanizawa1996

    <img src="schematic_float_Tanizawa1996.png" width="400px" />

    \cite{Tanizawa1997}と\cite{Tanizawa1996}には，BEM-MELを使う際の消波方法が紹介されている．
    この理論は簡単で，$`\frac{D\phi}{Dt}`$と$`\frac{D\bf x}{Dt}`$に減衰項を追加するだけである．
    興味深いのは，この項は，単なる消波だけでなく，ある領域の表面変位と速度ポテンシャルを理想的な値へと保つことができるため，造波装置側に設置するのも有効であることである．

    浮体の左右には整流板が設置されている．水深は4.5mで一定．

    | wave length | water height (m) | wave period from the linear dispersion relation (s) |
    |:-------|:------|:------|
    |  1.8 m  | 0.05 m | 1.07426 |
    |  2.7 m  | 0.05 m | 1.31570 |
    |  3.6 m  | 0.05 m | 1.51924 |
    |  1.8 m  | 0.1 m  | 1.07826 |    
    |  2.7 m  | 0.1 m  | 1.31570 |
    |  3.6 m  | 0.1 m  | 1.51924 |

    | Item of the floating body | M.K.S. system | Non-dimensional system |
    |:-------|:------|:-------|
    | Breadth | 0.74 m | 1.0 |
    | Draft | 0.25 m | 0.338 |
    | Displacement | 184.3 kg | 1.0 |
    | Center of gravity | 0.22 m | 0.297 |
    | Radius of inertia | 0.266 m | 0.359 |
    | Natural period of heave | 1.408 s | 5.12 |
    | Natural period of roll | 1.775 s | 6.46 |
    | Spring constant of mooning | 51.07 N/m | 0.00704 |

    '''

    start = 0.

    # L = 1.8
    a = H/2
    h = 4.5
    z_surface = h

    id = SimulationCase
    id += "_H" + str(H).replace(".", "d")

    if T is not None:
        id += "_T" + str(T).replace(".", "d")
    if L is not None:
        id += "_L" + str(L).replace(".", "d")

    if wavemaker_type == "":
        # wavemaker_type = "piston"
        wavemaker_type = "potential"
        # wavemaker_type = "flap"

    # if meshname == "":
    #     meshname = "water_no_float0d1"

    id = add_id(id)

    if "no_float" in suffix:
        h = 5
        z_surface = h

    water = {"name": "water", "type": "Fluid"}
    tank = {"name": "tank", "type": "RigidBody", "isFixed": True}
    absorber = {"name": "absorber", "type": "Absorber", "isFixed": True}

    if T is not None:
        absorber["wave_theory"] = [0, T, h, 0]
    if L is not None:
        absorber["wave_theory_L"] = [0, L, h, 0]

    gauges = []
    dx = 0.5
    for i in range(10):
        gauges.append({"name": "gauge"+str(i),
                       "type": "wave gauge",
                       "position": [-2.5 + dx*i, 0, h + 2, 0.2 + dx*i, 0, h - 2]})

    phase_shift = -pi/2.
    # if wavemaker_type contains "piston":

    if "piston" in wavemaker_type:
        wavemaker = {"name": "wavemaker",
                     "type": "RigidBody",
                     "velocity": ["piston_using_wave_length", start, a, L, h, 1, 0, 0]}
        if T is not None:
            wavemaker["velocity"] = ["piston", start, a, T, h, 1, 0, 0]
        elif L is not None:
            wavemaker["velocity"] = [
                "piston_using_wave_length", start, a, L, h, 1, 0, 0]
    elif "flap" in wavemaker_type:
        wavemaker = {"name": "wavemaker",
                     "type": "RigidBody",
                     "velocity": ["flap_using_wave_length", start, a, L, h, h, 0, 1, 0]}
        if T is not None:
            wavemaker["velocity"] = ["flap", start, a, T, h, h, 0, 1, 0]
        elif L is not None:
            wavemaker["velocity"] = [
                "flap_using_wave_length", start, a, L, h, h, 0, 1, 0]
    elif "potential" in wavemaker_type:
        wavemaker = {"name": "wavemaker",
                     "type": "SoftBody",
                     "isFixed": True,
                     "velocity": ["velocity_using_wave_length", start, a, L, h, h]}
        if T is not None:
            wavemaker["velocity"] = ["velocity", start, a, T, h, h]
        elif L is not None:
            wavemaker["velocity"] = [
                "velocity_using_wave_length", start, a, L, h, h]

    Lx = 0.74
    Lz = 0.415
    draft = 0.25
    Ly = 1

    floats = []
    if "multiple" in id:
        objfolder = code_home_dir + "/cpp/obj/Tanizawa1996_multiple"
        # i 1 to 8
        # j 1 to 3
        # COM = [12 - 1.5*(i-1),-2 + 2 *(j-1) ,z_surface - d + 0.22]
        for j in range(1, 4):
            for i in range(1, 9):
                float = {"name": "float"+str(i+8*(j-1)),
                         "type": "RigidBody",
                         "velocity": "floating",
                         "isFixed": [False, True, False, True, False, True]}
                m = 184.3
                float["mass"] = m
                Ixx = 1./12.*m*(Ly*Ly+Lz*Lz)
                Iyy = m*math.pow(0.266, 2)
                Izz = 1./12.*m*(Lx*Lx+Ly*Ly)
                float["COM"] = [12 - 1.5*(i-1), -2 + 2 * (j-1), z_surface - draft + 0.22]
                float["spring"] = [float["COM"][0], float["COM"][1], float["COM"][2], 51., 1000., 0.]
                float["MOI"] = [100.*Iyy, Iyy, 100.*Iyy]
                float["objfile"] = objfolder+"/float"+str(i+8*(j-1))+".obj"
                floats.append(float)
        water["objfile"] = objfolder + "/water_mod_coarse.obj"
    else:
        objfolder = code_home_dir + "/cpp/obj/Tanizawa1996"
        m = 184.3  # 184.3 kg
        Ixx = 1./12.*m*(Ly*Ly+Lz*Lz)
        K = 0.266
        Iyy = m * K * K  # 0.266 is the radius of inertia
        Izz = 1./12.*m*(Lx*Lx+Ly*Ly)
        float = {"name": "float",
                 "type": "RigidBody",
                 "velocity": "floating",
                 # "damping": [100, 100, 100, 0, 100, 0, 0., 1000000.],
                 "isFixed": [False, True, False, True, False, True],
                 "mass": m}

        absorber1 = {"name": "absorber1",
                     "type": "Absorber",
                     "isFixed": True,
                     "objfile": objfolder+"/absorber1.obj",
                     "wave_theory_L": [a, L, h, 0]}

        absorber2 = {"name": "absorber2",
                     "type": "Absorber",
                     "isFixed": True,
                     "objfile": objfolder+"/absorber2.obj",
                     "wave_theory": [0, L, h, 0]}

        water["objfile"] = objfolder + "/water.obj"

        # 密度1.18g/cm^3と質量から計算，{Ixx,Iyy,Izz}={0.3034752725,0.1967978007,0.3692731657}
        # 板厚0.08mmと質量から計算した結果，{Ixx,Iyy,Izz}={0.3320104413,0.2204711483,0.4018598521}

        float["COM"] = [0., 0., z_surface - draft + 0.22]
        anchor = [float["COM"][0] - 2, float["COM"][1], float["COM"][2]]
        kx = 51
        ky = 0.
        kz = 0.
        float["spring"] = [float["COM"][0], float["COM"][1], float["COM"][2], anchor[0],
                           # 実際は[51., 0., 0.]だがy方向にずれが生じるためこのようにした
                           anchor[1], anchor[2], kx, ky, kz]
        print("COM ", float["COM"], " mass ", float["COM"], " Ly ", Ly)
        float["MOI"] = [1000*Iyy, Iyy, 1000*Iyy]
        float["objfile"] = objfolder+"/float.obj"
        floats.append(float)

    # water["objfile"] = objfolder + "/water_mod_coarse.obj"

    if "no_float" in id:
        wavemaker["objfile"] = objfolder + "/wavemaker_no_float.obj"
        tank["objfile"] = objfolder + "/tank_no_float.obj"
        absorber["objfile"] = objfolder+"/absorber.obj"
        water["objfile"] = objfolder + "/" + meshname + ".obj"
        inputfiles = [tank, water, wavemaker, absorber]
    else:
        wavemaker["objfile"] = objfolder + "/wavemaker.obj"
        tank["objfile"] = objfolder + "/tank.obj"
        absorber["objfile"] = objfolder+"/absorber.obj"
        inputfiles = [tank, water, absorber1, absorber2]

    inputfiles += gauges
    if "no_float" not in id:
        inputfiles += floats

    setting = {"max_dt": dt,
               "end_time_step": 20000,
               "end_time": 30,
               "element": element,
               "ALE": ALE,
               "ALEPERIOD": ALEPERIOD}

    generate_input_files(inputfiles, setting, IO_dir, id)
elif "AkitaHarbor" in SimulationCase:

    r'''DOC_EXTRACT 2_1_0_AkitaHarbor
    '''

    start = 0.
    # L = 1.8
    a = 0.#H/2
    h = 10.
    T = 5.
    z_surface = h
    id = SimulationCase
    id += "_H" + str(H).replace(".", "d")
    if T is not None:
        id += "_T" + str(T).replace(".", "d")
    if L is not None:
        id += "_L" + str(L).replace(".", "d")
    wavemaker_type = "potential"
    id = add_id(id)
    water = {"name": "water", "type": "Fluid"}
    tank = {"name": "wall", "type": "RigidBody", "isFixed": True}
    absorber = {"name": "absorber",
                "type": "Absorber",
                "isFixed": True}

    gauges = []
    dx = 0.5
    wavemaker = {"name": "wavemaker",
                "type": "Absorber",
                "isFixed": True}

    bottom_in_z = 0
    a = 2
    wave_theory = [0, T, h, bottom_in_z, theta] if theta is not None else [0, T, h, bottom_in_z]
    absorber["wave_theory"] = wave_theory    
    wave_theory2 = [a, T, h, bottom_in_z, theta] if theta is not None else [a, T, h, bottom_in_z]
    wavemaker["wave_theory"] = wave_theory2

    objfolder = code_home_dir + "/cpp/obj/AkitaHarbor"
    water["objfile"] = objfolder + "/water.obj"
    wavemaker["objfile"] = objfolder + "/maker.obj"
    tank["objfile"] = objfolder + "/wall.obj"
    absorber["objfile"] = objfolder+"/absorber.obj"
    # inputfiles = [tank, water, wavemaker, absorber]
    inputfiles = [tank, water, wavemaker]

    setting = {"max_dt": dt,
               "end_time_step": 20000,
               "end_time": 300,
               "element": element,
               "ALE": ALE,
               "ALEPERIOD": ALEPERIOD}

    generate_input_files(inputfiles, setting, IO_dir, id)
elif "Cheng2018" in SimulationCase:

    start = 0.

    T = 1.2
    H = 0.04
    a = H/2  # Ren 2015 used H=[0.1(a=0.05), 0.03(a=0.06), 0.04(a=0.02)]
    h = 0.4

    if "meshA" in SimulationCase:
        id0 = "meshA"
    elif "meshB" in SimulationCase:
        id0 = "meshB"
    elif "meshC" in SimulationCase:
        id0 = "meshC"
    else:
        id0 = "meshA"
    # id0 = "mesh_D"
    # id0 = "_no"
    # id0 = "_multiple"

    wavemaker_type = "piston"
    # wavemaker_type = "potential"
    # wavemaker_type = "flap"

    id = SimulationCase

    id += "_H" + str(H).replace(".", "d")
    id += "_T" + str(T).replace(".", "d")
    id += "_" + wavemaker_type
    id += "_with_float"

    water = {"name": "water", "type": "Fluid"}
    tank = {"name": "tank", "type": "RigidBody", "isFixed": True}
    absorber = {"name": "absorber", "type": "Absorber", "isFixed": True}

    gauges = []
    dx = 0.5
    for i in range(10):
        gauges.append({"name": "gauge"+str(i),
                       "type": "wave gauge",
                       "position": [0.2 + dx*i, 0, 0.6, 0.2 + dx*i, 0, 0.1]})

    z_surface = 0.4
    phase_shift = -pi/2.
    # if wavemaker_type contains "piston":
    if "free_decay" in id:
        wavemaker = {"name": "wavemaker",
                     "type": "RigidBody",
                     "velocity": ["piston", start, 0, T, h, 1, 0, 0]}
    elif "piston" in wavemaker_type:
        wavemaker = {"name": "wavemaker",
                     "type": "RigidBody",
                     "velocity": ["piston", start, a, T, h, 1, 0, 0]}
    elif "flap" in wavemaker_type:
        wavemaker = {"name": "wavemaker",
                     "type": "RigidBody",
                     "velocity": ["flap", start, a, T, h, h, 0, 1, 0]}
    else:
        wavemaker = {"name": "wavemaker",
                     "type": "SoftBody",
                     "isFixed": True,
                     "velocity": ["velocity", start, a, T, h, z_surface, phase_shift]}

    Lx = 0.3
    Lz = 0.2
    d = 0.1
    # 500*Lx*Ly*Lz (浮体全質量) = 1000*Lx*Ly*d (排除される水の質量)
    Ly = 0.42

    float = {"name": "float",
             "type": "RigidBody",
             "velocity": "floating",
             "isFixed": [False, True, False, True, False, True]}

    float["mass"] = m = 500*Lx*Lz*Ly
    Ixx = 1./12.*m*(Ly*Ly+Lz*Lz)
    # Iyy = 1./12.*m*(Lx*Lx+Lz*Lz)
    # 密度1.18g/cm^3と質量から計算，{Ixx,Iyy,Izz}={0.3034752725,0.1967978007,0.3692731657}
    Iyy = 0.1967978007
    # Iyy = 0.2204711483 # 板厚0.08mmと質量から計算した結果，{Ixx,Iyy,Izz}={0.3320104413,0.2204711483,0.4018598521}
    # Iyy = 0.15
    Izz = 1./12.*m*(Lx*Lx+Ly*Ly)
    z_surface = 0.4

    # 密度1.18g/cm^3と質量から計算，{Ixx,Iyy,Izz}={0.3034752725,0.1967978007,0.3692731657}
    # 板厚0.08mmと質量から計算した結果，{Ixx,Iyy,Izz}={0.3320104413,0.2204711483,0.4018598521}

    float["COM"] = [3.35, 0., z_surface - d]
    print("COM ", float["COM"], " mass ", float["COM"], " Ly ", Ly)
    float["MOI"] = [10**10*Ixx, Iyy, 10**10*Izz]

    if "free_decay" in id:
        objfolder = code_home_dir + "/cpp/obj/Cheng2018"
        water["objfile"] = objfolder + "/water_inclined_mod.obj"
        inputfiles = [tank, wavemaker, water, float, absorber]
        float["objfile"] = objfolder+"/float_inclined.obj"
        wavemaker["velocity"] = ["velocity", start,
                                 0., T, h, z_surface, phase_shift]
    elif "without_float" in id:
        objfolder = code_home_dir + "/cpp/obj/Cheng2018"
        water["objfile"] = objfolder + "/water_without_float9_mod.obj"
        inputfiles = [tank, wavemaker, water, absorber]
        float["objfile"] = objfolder+"/float7.obj"
    else:
        objfolder = code_home_dir + "/cpp/obj/Cheng2018"
        if id0 != "":
            water["objfile"] = objfolder + "/water8_" + id0 + ".obj"
        else:
            water["objfile"] = objfolder + "/water8_modmodmod.obj"
        inputfiles = [tank, wavemaker, water, float, absorber]
        float["objfile"] = objfolder+"/float7.obj"

    wavemaker["objfile"] = objfolder + "/wavemaker8.obj"
    tank["objfile"] = objfolder + "/tank4.obj"

    absorber["objfile"] = objfolder+"/absorber.obj"
    inputfiles += gauges

    setting = {"max_dt": 0.02,
               "end_time_step": 10000,
               "end_time": 20}

    generate_input_files(inputfiles, setting, IO_dir, id)
elif "Ren2015" in SimulationCase:

    r'''DOC_EXTRACT 2_1_0_validation_Ren2015

    <img src="schematic_Ren2015.png" width="400px" />

    This case based on \cite{Ren2015} is for the validation of the floating body motion analysis using the BEM-MEL.
    The floating body is a rectangular box with the dimension of $`(l_x,l_y,l_z)=(0.3,0.42,0.2) {\rm m}`$ 
    The density of the floating body is $`0.5\times1000 {\rm kg/m^3}`$.
    The moment of inertia of the floating body is $`(I_{xx},I_{yy},I_{zz}) = (\frac{m}{12}(l_y^2+l_z^2),\frac{m}{12}(l_x^2+l_z^2),\frac{m}{12}(l_x^2+l_y^2))`$.

    You can find numerical results compared with this case from \cite{Cheng2018} and \cite{Bihs2017}.

    [Youtube DualSPHysics](https://www.youtube.com/watch?v=VDa4zcMDjJA)

    '''

    start = 0.

    T = 1.2
    if H is None:
        H = 0.04
    a = H/2  # Ren 2015 used H=[0.1(a=0.05), 0.03(a=0.06), 0.04(a=0.02)]
    h = 0.4

    id0 = ""
    # id0 = "_no"
    # id0 = "_multiple"

    wavemaker_type = "piston"
    # wavemaker_type = "potential"

    id = SimulationCase + id0 + "_H"+str(H).replace(".", "d")
    id += "_T"+str(T).replace(".", "d")

    if theta is not None:
        id += "_theta"+str(theta).replace(".", "d")

    water = {"name": "water", "type": "Fluid"}
    tank = {"name": "tank",
            "type": "RigidBody",
            "isFixed": True}

    z_surface = 0.4

    if wavemaker_type == "piston":
        wavemaker = {"name": "wavemaker",
                     "type": "RigidBody",
                     "velocity": ["piston", start, a, T, h, 1, 0, 0]}
    else:
        wavemaker = {"name": "wavemaker",
                     "type": "SoftBody",
                     "isFixed": True,
                     "velocity": ["linear_traveling_wave", start, a, T, h, z_surface]}

    L = 0.3
    d = 0.1
    # 500*L*W*H (浮体全質量) = 1000*L*W*d (排除される水の質量)
    # d = H/2 = 0.1
    W = 0.5

    float = {"name": "float",
             "type": "RigidBody",
             "velocity": "floating"}

    float_H = 0.2
    float_W = 0.42
    float["mass"] = m = 500*L*float_H*float_W
    Ixx = 1./12.*m*(float_W*float_W+float_H*float_H)
    Iyy = 1./12.*m*(L*L+float_H*float_H)
    # Iyy = 0.13
    Izz = 1./12.*m*(L*L+float_W*float_W)
    MOI = [Ixx, Iyy, Izz]
    z_surface = 0.4
    z_floatinbody_bottom = z_surface - d
    # float["COM"] = [2+L/2., W/2, z_surface]
    # float["COM"] = [4.6+L/2, W/2, z_surface]
    float["COM"] = [0, W/2, z_surface]
    print("COM ", float["COM"], " mass ", float["COM"], " W ", W)
    float["MOI"] = [10**10*Ixx, Iyy, 10**10*Izz]

    # if id contains "multiple":
    if "multiple" in id:
        float["COM"] = [4.3+L/2, W/2, z_surface]
        objfolder = code_home_dir + "/cpp/obj/Ren2015_multiple"

        # Initialize the object files
        water["objfile"] = f"{objfolder}/water400.obj"
        wavemaker["objfile"] = f"{objfolder}/wavemaker20.obj"
        tank["objfile"] = f"{objfolder}/tank.obj"
        float["objfile"] = f"{objfolder}/float10.obj"

        inputfiles = [tank, wavemaker, water]

        # Create float_a to float_g
        float_names = ['a', 'b', 'c', 'd', 'e', 'f', 'g']
        float_offsets = [-3, -2, -1, 0, 1, 2, 3]

        for name, offset in zip(float_names, float_offsets):
            new_float = copy.deepcopy(float)
            new_float["name"] = f"float_{name}"
            new_float["objfile"] = f"{objfolder}/float_{name}50.obj"
            new_float["COM"][0] = float["COM"][0] + offset
            inputfiles.append(new_float)
    elif "_no" in id:
        objfolder = code_home_dir + "/cpp/obj/Ren2015_no_float"
        water["objfile"] = objfolder + "/water25.obj"
        wavemaker["objfile"] = objfolder + "/wavemaker20.obj"
        tank["objfile"] = objfolder + "/tank2.obj"
        inputfiles = [tank, wavemaker, water]
    elif "2D" in id:
        float["COM"] = [0, 0, z_surface]  # 新しい設定
        objfolder = code_home_dir + "/cpp/obj/Ren2015_2D"
        water["objfile"] = objfolder + "/water.obj"
        wavemaker["objfile"] = objfolder + "/wavemaker.obj"
        tank["objfile"] = objfolder + "/tank.obj"
        float["objfile"] = objfolder+"/float.obj"
        bottom_in_z = 0
        a = H / 2
        h = 0.4
        absorber1 = {"name": "absorber1",
                     "type": "Absorber",
                     "isFixed": True,
                     "objfile": objfolder+"/absorber1.obj",
                     "wave_theory": [a, T, h, bottom_in_z]}

        absorber2 = {"name": "absorber2",
                     "type": "Absorber",
                     "isFixed": True,
                     "objfile": objfolder+"/absorber2.obj",
                     "wave_theory": [0, T, h, bottom_in_z]}

        # inputfiles = [tank, wavemaker, water, float, absorber]
        inputfiles = [tank, water, float, absorber1, absorber2]
    else:
        # float["COM"] = [5, 0, z_surface]#新しい設定
        float["COM"] = [0, 0, z_surface]  # 新しい設定
        objfolder = code_home_dir + "/cpp/obj/Ren2015"
        # water["objfile"] = objfolder + "/water400meshlab.obj"
        water["objfile"] = objfolder + "/water.obj"
        wavemaker["objfile"] = objfolder + "/wavemaker.obj"
        tank["objfile"] = objfolder + "/tank.obj"
        float["objfile"] = objfolder+"/float.obj"
        bottom_in_z = 0
        a = H / 2
        h = 0.4
        wave_theory = [a, T, h, bottom_in_z,
                       theta] if theta is not None else [a, T, h, bottom_in_z]
        absorber = {"name": "absorber",
                    "type": "Absorber",
                    "isFixed": True,
                    "objfile": objfolder+"/absorber.obj",
                    "wave_theory": wave_theory}
        # inputfiles = [tank, wavemaker, water, float, absorber]
        inputfiles = [tank, water, float, absorber]

    id = add_id(id)

    setting = {"max_dt": dt,
               "end_time_step": 10000,
               "end_time": 10,
               "element": element,
               "ALE": ALE,
               "ALEPERIOD": ALEPERIOD}

    generate_input_files(inputfiles, setting, IO_dir, id)
elif "Hadzic2005" in SimulationCase:

    r'''DOC_EXTRACT 2_1_0_validation_Hadzic2005                

    ## Hadzic2005

    <img src="schematic_Hadzic2005.png" width="400px"/>

    This case based on \cite{Hadzic2005} is for the validation of the floating body motion analysis using the BEM-MEL.        
    The floating body is a rectangular box with the dimension of L10 cm x H5 cm x W29 cm.        
    The density of the floating body is 0.68x1000 kg/m^3, therefore the mass of the floating body is 0.68x0.05x0.1x0.29x1000 kg.
    The moment of inertia of the floating body is 14 kg cm^2.

    [CAD data](https://a360.co/46CisV7)

    [spheric Test 12](https://www.spheric-sph.org/tests/test-12)

    [Youtube Nextflow](https://www.youtube.com/watch?v=H92xupH9508)

    '''

    # id = SimulationCase + "_multiple"
    # id = SimulationCase + "_without_float"
    id = SimulationCase

    if "multiple" in id:
        objfolder = code_home_dir + "/cpp/obj/Hadzic2005_24floats"
    elif "without_float" in id:
        objfolder = code_home_dir + "/cpp/obj/Hadzic2005"
    else:
        objfolder = code_home_dir + "/cpp/obj/Hadzic2005"

    water = {"name": "water",
             "type": "Fluid"}

    if "multiple" in id:
        water["objfile"] = objfolder + "/water1000meshlab.obj"
    elif "without_float" in id:
        water["objfile"] = objfolder + "/water_without_float9.obj"
    else:
        water["objfile"] = objfolder + "/water1000.obj"

    vertices, triangles = read_obj(water["objfile"])
    water["vertices"] = vertices
    water["triangles"] = triangles

    tank = {"name": "tank",
            "type": "RigidBody",
            "isFixed": True,
            "objfile": objfolder + "/tank10.obj"}

    start_time = 0.

    wavemaker = {"name": "wavemaker",
                 "type": "RigidBody",
                 "objfile": objfolder + "/wavemaker100.obj",
                 "velocity": ["Hadzic2005", start_time],
                 "acceleration": ["Hadzic2005", start_time],
                 "COM": [0., 0., 0.]}

    float = {"name": "float",
             "type": "RigidBody",
             "objfile": objfolder+"/float20.obj",
             "output": "json",
             "velocity": "floating"}

    L = 0.1
    W = 0.29
    H = 0.05
    A = L*W
    # d = 0.03 #喫水深さ draft depth (draught)
    density = 680
    float["mass"] = m = density*0.05*0.1*0.29
    d = m / (rho * A)
    print("d", d)
    # MOI = 14.*(0.01*0.01)  # original kg*m*m
    Ixx = 1./12.*m*(W*W+H*H)
    Iyy = 1./12.*m*(L*L+H*H)
    Izz = 1./12.*m*(L*L+W*W)
    z_surface = 0.4
    z_floatinbody_bottom = z_surface - d
    # float["COM"] = [-(4.-2.11), 0., z_floatinbody_bottom + 0.05/2]
    float["COM"] = [2.11, W/2, z_floatinbody_bottom + H/2]
    # float["MOI"] = [Ixx*10**10,Iyy,Izz*10**10]
    float["MOI"] = [Ixx*10**10, 0.001589233395, Izz*10**10]
    # float["MOI"] = [Ixx*10**10,0.0019,Izz*10**10]

    if "without_float" in id:
        inputfiles = [tank, wavemaker, water]
    else:
        inputfiles = [tank, wavemaker, water, float]

    setting = {"max_dt": 0.05,
               "end_time_step": 10000,
               "end_time": 9}

    generate_input_files(inputfiles, setting, IO_dir, id)
elif "Kramer2021" in SimulationCase:

    r'''DOC_EXTRACT 2_1_1_validation_Kramer2021

    ## Kramer2021

    This case is for the validation of the floating body motion analysis using the BEM-MEL.

    <img src="schematic_Kramer2021.png" width="400px" />

    The floating body is a sphere with the diameter of 0.3 m.
    The mass of the floating body is 7.056 kg.
    The moment of inertia of the floating body is set to be almost infinite to ignore the effect of the rotation.

    The sphere is dropped from the height of 0.03 m above the water surface.
    '''

    D = 300/1000
    # id = "H0"+str(H0).replace(".", "d")

    water = {"name": "water", "type": "Fluid"}
    tank = {"name": "tank", "type": "RigidBody", "isFixed": True}

    start_t = 0.02
    float = {"name": "float",
             "type": "RigidBody",
             "velocity": ["floating", start_t]}

    float["mass"] = m = 7.056
    # float["reverseNormal"] = True
    z_surface = 900/1000
    float["radius_of_gyration"] = [10**10, 10**10, 10**10]
    float["MOI"] = [m*math.pow(float["radius_of_gyration"][0], 2),
                    m*math.pow(float["radius_of_gyration"][1], 2),
                    m*math.pow(float["radius_of_gyration"][2], 2)]
    # float["translate"] = [0., 0., 0.1*D + 900/1000]
    z_surface = 900/1000
    id = SimulationCase
    if "H00d03" in id or "H00d03" in suffix:
        H0 = D*0.1
        float["COM"] = [0., 0., z_surface + H0]  # 今回は重要ではない
        objfolder = code_home_dir + "/cpp/obj/Kramer2021_H00d03"
        water["objfile"] = objfolder + "/water.obj"
        tank["objfile"] = objfolder + "/tank.obj"
        float["objfile"] = objfolder + "/sphere.obj"
    elif "H00d09" in id or "H00d09" in suffix:
        H0 = D*0.3
        float["COM"] = [0., 0., z_surface + H0]  # 今回は重要ではない
        objfolder = code_home_dir + "/cpp/obj/Kramer2021_H00d09"
        water["objfile"] = objfolder + "/water.obj"
        tank["objfile"] = objfolder + "/tank.obj"
        float["objfile"] = objfolder + "/sphere.obj"
    else:
        # throw error
        print("Error: SimulationCase is not defined")
        exit()

    id = add_id(id)

    gauges = []
    dx = 0.6
    for i in range(10):
        gauges.append({"name": "gauge"+str(i), "type": "wave gauge",
                      "position": [0. + dx*i, 0, 1.1, 0. + dx*i, 0, 1.1-0.6]})

    inputfiles = [tank, water, float]
    inputfiles += gauges

    rho = 998.2
    g = 9.82

    setting = {"max_dt": dt,
               "WATER_DENSITY": rho,
               "GRAVITY": g,
               "end_time_step": 1000000,
               "end_time": 15,
               "element": element,
               "ALE": ALE,
               "ALEPERIOD": ALEPERIOD}

    generate_input_files(inputfiles, setting, IO_dir, id)
elif "Palm2016" in SimulationCase:

    r'''DOC_EXTRACT 2_1_0_validation_Palm2016

    ## Palm2016

    | wave height (m) | wave period (s) |
    |:-------|:------|
    | 0.04   | 1.0   |
    | 0.04   | 1.2   |
    | 0.04   | 1.4   |
    | 0.08   | 1.0   |
    | 0.08   | 1.2   |
    | 0.08   | 1.4   |

    '''

    # % 造波機の設定
    start = 0.

    if H is None:
        H = 0.04

    a = H/2
    T = 1.
    h = 0.9
    z_surface = 0.9

    # SimulationCase += "_H" + str(H).replace(".", "d")
    # SimulationCase += "_T" + str(T).replace(".", "d")

    id = add_id(id)

    rho = 1000
    draft = 0.172

    # % 浮体の設定
    D = 0.515  # m (直径)
    vol = math.pi*D**2./4. * draft
    M = rho*vol  # kg (質量)
    Ixx = 0.9  # kg*m^2 (慣性モーメント)
    float_bottom = z_surface - draft  # m (底面の高さ)
    COM = [6, 0, 0.0758 + float_bottom]  # m (重心位置)

    # % 係留索の設定  120度の角度で3本の係留索を設定
    i = 0
    horizontal_length = 1.66
    r = D/2 + 0.015
    q = i*2*pi/3
    X_fair_leadA = [r * math.cos(q) + COM[0], r *
                    math.sin(q) + COM[1], z_surface]
    X_anchorA = Add3d(X_fair_leadA, [
                      horizontal_length*math.cos(q), horizontal_length*math.sin(q), -z_surface])

    q = (i+1)*2*pi/3
    X_fair_leadB = [r * math.cos(q) + COM[0], r *
                    math.sin(q) + COM[1], z_surface]
    X_anchorB = Add3d(X_fair_leadB, [
                      horizontal_length*math.cos(q), horizontal_length*math.sin(q), -z_surface])

    q = (i+2)*2*pi/3
    X_fair_leadC = [r * math.cos(q) + COM[0], r *
                    math.sin(q) + COM[1], z_surface]
    X_anchorC = Add3d(X_fair_leadC, [
                      horizontal_length*math.cos(q), horizontal_length*math.sin(q), -z_surface])

    total_lengthA = Norm3d(Differece3d(X_anchorA, X_fair_leadA))
    total_lengthB = Norm3d(Differece3d(X_anchorB, X_fair_leadB))
    total_lengthC = Norm3d(Differece3d(X_anchorC, X_fair_leadC))

    stiffness = 300*10**6  # ! [N/m]
    damp = .5  # ! [N/(m/s^2)]
    density = 0.1447  # ! [kg/m]
    n_points = 30
    diam = 4.786/1000
    stiffness = stiffness * (math.pi*(diam/2.)**2)
    water = {"name": "water",
             "type": "Fluid"}

    tank = {"name": "tank",
            "type": "RigidBody",
            "isFixed": True}

    wavemaker = {"name": "wavemaker",
                 "type": "SoftBody",
                 "isFixed": True,
                 "velocity": ["velocity", start, a, T, h, z_surface]}

    float = {"name": "float",
             "type": "RigidBody",
             "velocity": "floating",
             "mass": M,
             "COM": COM,
             "MOI": [Ixx, Ixx, Ixx*1000.]}

    id = SimulationCase

    # if meshname == "":
    #     meshname = "water_mod"
    # id += "_MESH" + meshname

    # if dt is not None:
    #     max_dt = dt
    # else:
    #     max_dt = 1/20

    # id += "_DT" + str(max_dt).replace(".", "d")

    # if element != "":
    #     id += "_ELEM" + element

    # if ALE != "":
    #     id += "_ALE" + ALE

    # if ALEPERIOD != "":
    #     id += "_ALEPERIOD" + ALEPERIOD

    # if suffix != "":
    #     id += "_" + suffix

    if "with_mooring" in id:
        float["mooringA"] = ["mooringA",
                             X_anchorA[0], X_anchorA[1], X_anchorA[2],
                             X_fair_leadA[0], X_fair_leadA[1], X_fair_leadA[2],
                             total_lengthA,
                             n_points,
                             density,
                             stiffness,
                             damp,
                             diam]

        float["mooringB"] = ["mooringB",
                             X_anchorB[0], X_anchorB[1], X_anchorB[2],
                             X_fair_leadB[0], X_fair_leadB[1], X_fair_leadB[2],
                             total_lengthB,
                             n_points,
                             density,
                             stiffness,
                             damp,
                             diam]

        float["mooringC"] = ["mooringC",
                             X_anchorC[0], X_anchorC[1], X_anchorC[2],
                             X_fair_leadC[0], X_fair_leadC[1], X_fair_leadC[2],
                             total_lengthC,
                             n_points,
                             density,
                             stiffness,
                             damp,
                             diam]

    probe1 = {"name": "probe1",
              "type": "wave gauge",
              "position": [float["COM"][0] - 0.6, float["COM"][1] + 2., 1.2, float["COM"][0] - 0.6, float["COM"][1] + 2., 0.6]}
    probe2 = {"name": "probe2",
              "type": "wave gauge",
              "position": [float["COM"][0], float["COM"][1] + 2., 1.2, float["COM"][0], float["COM"][1] + 2., 0.6]}
    probe3 = {"name": "probe3",
              "type": "wave gauge",
              "position": [float["COM"][0] + 0.14, float["COM"][1] + 2., 1.2, float["COM"][0] - 0.14, float["COM"][1] + 2., 0.6]}
    probe4 = {"name": "probe4",
              "type": "wave gauge",
              "position": [float["COM"][0] - 0.31, float["COM"][1] + 2., 1.2, float["COM"][0] - 0.31, float["COM"][1] + 2., 0.6]}

    probes = [probe1, probe2, probe3, probe4]

    objfolder = code_home_dir + "/cpp/obj/Palm2016"
    # water["objfile"] = objfolder + "/water_mod.obj"
    water["objfile"] = objfolder + "/" + meshname + ".obj"
    wavemaker["objfile"] = objfolder + "/wavemaker_mod.obj"
    tank["objfile"] = objfolder + "/tank_mod.obj"
    float["objfile"] = objfolder+"/float_mod.obj"

    inputfiles = [tank, wavemaker, water, float]
    inputfiles += probes

    setting = {"max_dt": dt,
               "end_time_step": 100000,
               "end_time": 30,
               "element": element,
               "ALE": ALE,
               "ALEPERIOD": ALEPERIOD}

    generate_input_files(inputfiles, setting, IO_dir, id)
elif "Liang2022" in SimulationCase:

    r'''DOC_EXTRACT 2_1_0_validation_Liang2022

    ## Liang2022

    \cite{Liang2022}
    Shandong Provincial Key Laboratory of Ocean Engineering, Ocean University of China.
    The flume was 3 m wide, 60 m long, and 1.8 m deep.
    The flume is divide into two, 0.8 m and 2.2 m wide, and this test was conducted in the 0.8 m wide flume.

    wave period range: 1 - 1.8 s
    wave height range: 0.05 - 0.15 m

    | wave height (m) | wave period (s) | water depth (m) |
    |:-------|:------|:-------|
    |  0.05  | 1.2   | 0.6    |
    |  0.07  | 1, 1.1, 1.2 ,1.4, 1.6, 1.8 |     |
    |  0.10  | 1, 1.1, 1.2 ,1.4, 1.6, 1.8 |     |
    |  0.13  | 1, 1.1, 1.2 ,1.4, 1.6, 1.8 |     |
    |  0.15  | 1.2   |     |

    The floating breakwater:

    | width | length | height | draft | mass | MOI | COM |
    |:---|:---|:---|:---|:---|:---|:---|
    | 0.5 m | 0.745 m | 0.28 m | 0.16 m | 58.09 kg | 2.441 kg m^2 | 0.0652 m from the bottom |

    | types | chain length | stiffness |
    |:---|:---|:---|
    | type A | 1.567 m | 2.36*10^3 N/m |
    | type B | 1.125 m | 2.36*10^3 N/m |
    | type C | 0.809 m | 2.36*10^3 N/m |
    | type D | None | Fixed |

    The floating breakwater was made of 15 mm thick acrylic plates.

    Mooring system:
    The mooring line was made of  stainless steel with a line density of 0.177 kg/m.

    The wave gauges were WG1: 3.5 m from the front of the float, WG2: 3.0 m from the front of the float, 
    WG3: 3.0 m from the rear of the float, and WG4: 3.5 m from the rear of the float.

    '''

    # % 造波機の設定
    start = 0.
    H = 0.07
    a = H/2
    T = 1.4
    h = 0.6
    z_surface = h

    SimulationCase += "_H" + str(H).replace(".", "d")
    SimulationCase += "_T" + str(T).replace(".", "d")

    rho = 1000
    draft = 0.16
    Lx = W = 0.5  # float width in x direction
    Ly = 0.745

    V = draft * Lx * Ly
    print("V ", V, ", mass should be ", V*rho)

    # % 浮体の設定
    M = 58.09  # kg (質量)
    Iyy = 2.441  # kg*m^2 (慣性モーメント)
    float_bottom = z_surface - draft  # m (底面の高さ)
    COM = [4.5+0.5/2, 0, 0.0652 + float_bottom]  # m (重心位置)

    water = {"name": "water",
             "type": "Fluid"}

    tank = {"name": "tank",
            "type": "RigidBody",
            "isFixed": True}

    wavemaker = {"name": "wavemaker",
                 "type": "RigidBody",
                 "velocity": ["piston", start, a, T, h, 1, 0, 0]}

    float = {"name": "float",
             "type": "RigidBody",
             "velocity": "floating",
             "mass": M,
             "COM": COM,
             "MOI": [10.**10, Iyy, 10.**10]}

    # % 係留索の設定
    i = 0

    if "typeA" in SimulationCase:
        r = math.sqrt(1.567**2 - float_bottom**2)  # horizontal_length
        X_fair_leadA = [float["COM"][0] - Lx/2,
                        float["COM"][1] - Ly/2, float_bottom]
        X_anchorA = [float["COM"][0] - Lx/2 + r, float["COM"][1] - Ly/2, 0]
        X_fair_leadB = [float["COM"][0] - Lx/2,
                        float["COM"][1] + Ly/2, float_bottom]
        X_anchorB = [float["COM"][0] - Lx/2 + r, float["COM"][1] + Ly/2, 0]
        X_fair_leadC = [float["COM"][0] + Lx/2,
                        float["COM"][1] - Ly/2, float_bottom]
        X_anchorC = [float["COM"][0] + Lx/2 - r, float["COM"][1] - Ly/2, 0]
        X_fair_leadD = [float["COM"][0] + Lx/2,
                        float["COM"][1] + Ly/2, float_bottom]
        X_anchorD = [float["COM"][0] + Lx/2 - r, float["COM"][1] + Ly/2, 0]
    elif "typeB" in SimulationCase:
        r = math.sqrt(1.125**2 - float_bottom**2)
        X_fair_leadA = [float["COM"][0] - Lx/2,
                        float["COM"][1] - Ly/2, float_bottom]
        X_anchorA = [float["COM"][0] - Lx/2 - r, float["COM"][1] - Ly/2, 0]
        X_fair_leadB = [float["COM"][0] - Lx/2,
                        float["COM"][1] + Ly/2, float_bottom]
        X_anchorB = [float["COM"][0] - Lx/2 - r, float["COM"][1] + Ly/2, 0]
        X_fair_leadC = [float["COM"][0] + Lx/2,
                        float["COM"][1] - Ly/2, float_bottom]
        X_anchorC = [float["COM"][0] + Lx/2 + r, float["COM"][1] - Ly/2, 0]
        X_fair_leadD = [float["COM"][0] + Lx/2,
                        float["COM"][1] + Ly/2, float_bottom]
        X_anchorD = [float["COM"][0] + Lx/2 + r, float["COM"][1] + Ly/2, 0]
    elif "typeC" in SimulationCase:
        r = math.sqrt(0.809**2 - float_bottom**2)
        X_fair_leadA = [float["COM"][0] - Lx/2,
                        float["COM"][1] - Ly/2, float_bottom]
        X_anchorA = [float["COM"][0] - Lx/2 - r, float["COM"][1] - Ly/2, 0]
        X_fair_leadB = [float["COM"][0] - Lx/2,
                        float["COM"][1] + Ly/2, float_bottom]
        X_anchorB = [float["COM"][0] - Lx/2 - r, float["COM"][1] + Ly/2, 0]
        X_fair_leadC = [float["COM"][0] + Lx/2,
                        float["COM"][1] - Ly/2, float_bottom]
        X_anchorC = [float["COM"][0] + Lx/2 + r, float["COM"][1] - Ly/2, 0]
        X_fair_leadD = [float["COM"][0] + Lx/2,
                        float["COM"][1] + Ly/2, float_bottom]
        X_anchorD = [float["COM"][0] + Lx/2 + r, float["COM"][1] + Ly/2, 0]

    stiffness = 2.36*10**3  # ! [N/m]
    damp = .4  # ! [N/(m/s^2)]
    density = 0.177  # ! [kg/m]

    total_lengthA = Norm3d(Differece3d(X_anchorA, X_fair_leadA))
    total_lengthB = Norm3d(Differece3d(X_anchorB, X_fair_leadB))
    total_lengthC = Norm3d(Differece3d(X_anchorC, X_fair_leadC))
    total_lengthD = Norm3d(Differece3d(X_anchorD, X_fair_leadD))

    n_points = 30
    diam = 1*0.1  # 1 cm

    float["mooringA"] = ["mooringA",
                         X_anchorA[0], X_anchorA[1], X_anchorA[2],
                         X_fair_leadA[0], X_fair_leadA[1], X_fair_leadA[2],
                         total_lengthA,
                         n_points,
                         density,
                         stiffness,
                         damp,
                         diam]

    float["mooringB"] = ["mooringB",
                         X_anchorB[0], X_anchorB[1], X_anchorB[2],
                         X_fair_leadB[0], X_fair_leadB[1], X_fair_leadB[2],
                         total_lengthB,
                         n_points,
                         density,
                         stiffness,
                         damp,
                         diam]

    float["mooringC"] = ["mooringC",
                         X_anchorC[0], X_anchorC[1], X_anchorC[2],
                         X_fair_leadC[0], X_fair_leadC[1], X_fair_leadC[2],
                         total_lengthC,
                         n_points,
                         density,
                         stiffness,
                         damp,
                         diam]

    float["mooringD"] = ["mooringD",
                         X_anchorD[0], X_anchorD[1], X_anchorD[2],
                         X_fair_leadD[0], X_fair_leadD[1], X_fair_leadD[2],
                         total_lengthD,
                         n_points,
                         density,
                         stiffness,
                         damp,
                         diam]

    id = SimulationCase

    probe1 = {"name": "WG1",
              "type": "wave gauge",
              "position": [float["COM"][0] - W/2 - 3.5, float["COM"][1], 1.2, float["COM"][0] - W/2 - 3.5, float["COM"][1], 0.3]}
    probe2 = {"name": "WG2",
              "type": "wave gauge",
              "position": [float["COM"][0] - W/2 - 3.0, float["COM"][1], 1.2, float["COM"][0] - W/2 - 3.0, float["COM"][1], 0.3]}
    probe3 = {"name": "WG3",
              "type": "wave gauge",
              "position": [float["COM"][0] + W/2 + 3.0, float["COM"][1], 1.2, float["COM"][0] + W/2 + 3.0, float["COM"][1], 0.3]}
    probe4 = {"name": "WG4",
              "type": "wave gauge",
              "position": [float["COM"][0] + W/2 + 3.5, float["COM"][1], 1.2, float["COM"][0] + W/2 + 3.5, float["COM"][1], 0.3]}
    probes = [probe1, probe2, probe3, probe4]

    objfolder = code_home_dir + "/cpp/obj/Liang2022"
    water["objfile"] = objfolder + "/water5_mod.obj"
    wavemaker["objfile"] = objfolder + "/wavemaker5.obj"
    tank["objfile"] = objfolder + "/tank5.obj"
    float["objfile"] = objfolder+"/float5.obj"

    inputfiles = [tank, wavemaker, water, float]
    inputfiles += probes

    setting = {"max_dt": 0.03,
               "end_time_step": 100000,
               "end_time": 30}

    input_dir += SimulationCase
    os.makedirs(input_dir, exist_ok=True)
    output_dir = sys_home_dir + "/BEM/" + SimulationCase
    os.makedirs(output_dir, exist_ok=True)
    generate_input_files(inputfiles, setting, IO_dir, id)
elif "WaveGeneration" in SimulationCase:
    # objfolder = code_home_dir + "/cpp/obj/WaveGeneration"

    objfolder = code_home_dir + "/cpp/obj/Hadzic2005"

    water = {"name": "water",
             "type": "Fluid",
             "objfile": objfolder + "/water_without_float9_mod.obj"}

    vertices, triangles = read_obj(water["objfile"])
    water["vertices"] = vertices
    water["triangles"] = triangles

    tank = {"name": "tank", "type":
            "RigidBody", "isFixed": True,
            "objfile": objfolder + "/tank10.obj"}

    start = 0.

    T = 1.
    H = 0.05
    a = H/2  # Ren 2015 used H=[0.1(a=0.05), 0.03(a=0.06), 0.04(a=0.02)]
    h = 0.4

    # wavemaker_type = "piston"
    # wavemaker_type = "flap"
    # wavemaker_type = "potential"

    id = SimulationCase
    id += "_H" + str(H).replace(".", "d")
    id += "_T"+str(T).replace(".", "d")
    id += "_h"+str(h).replace(".", "d")

    z_surface = h

    if "piston" in SimulationCase:
        wavemaker = {"name": "wavemaker",
                     "type": "RigidBody",
                     "velocity": ["piston", start, a, T, h, 1, 0, 0],
                     "objfile": f"{objfolder}/wavemaker10.obj"}

    elif "flap" in SimulationCase:
        wavemaker = {"name": "wavemaker",
                     "type": "RigidBody",
                     "velocity": ["flap", start, a, T, h, h, 0, 1, 0],
                     "objfile": f"{objfolder}/wavemaker10.obj"}
    else:
        wavemaker = {"name": "wavemaker",
                     "type": "SoftBody",
                     "isFixed": True,
                     "velocity": ["velocity", start, a, T, h, z_surface],
                     "objfile": f"{objfolder}/wavemaker10.obj"}

    inputfiles = [tank, wavemaker, water]

    setting = {"max_dt": 0.03,
               "end_time_step": 10000,
               "end_time": 9}

    generate_input_files(inputfiles, setting, IO_dir, id)
elif "fish_without_free_surface" in SimulationCase:
    start = 0.

    id = SimulationCase

    objfolder = code_home_dir + "/cpp/obj/fish"

    waterA = {"name": "waterA",
              "type": "Fluid",
              "reverseNormal": True,
              "objfile": objfolder + "/bodyA20.obj"}

    waterB = {"name": "waterB",
              "type": "Fluid",
              "reverseNormal": True,
              "objfile": objfolder + "/bodyB20.obj"}

    waterC = {"name": "waterC",
              "type": "Fluid",
              "reverseNormal": True,
              "objfile": objfolder + "/bodyC20.obj"}

    bodyA = {"name": "bodyA",
             "type": "RigidBody",
             "COM": [0., 0, 0.25],
             "mass": 10**10,
             "MOI": [10**10, 10**10, 10**10],
             "output": "json",
             #  "velocity": ["sin", 0, 0.1, 5, 0, 0, 0, 0, 0, 1],
             "velocity": ["file", "./study_fish/bodyA.dat"],
             "objfile": objfolder + "/bodyA20.obj"}

    bodyB = {"name": "bodyB",
             "type": "RigidBody",
             "COM": [0.35, 0, 0.25],
             "mass": 10**10,
             "MOI": [10**10, 10**10, 10**10],
             "output": "json",
             #  "velocity": ["sin", 0, 0.1, 5, 0, 0, 0, 0, 0, 1],
             "velocity": ["file", "./study_fish/bodyB.dat"],
             "objfile": objfolder + "/bodyB20.obj"}

    bodyC = {"name": "bodyC",
             "type": "RigidBody",
             "COM": [0.7, 0, 0.25],
             "mass": 10**10,
             "MOI": [10**10, 10**10, 10**10],
             "output": "json",
             #  "velocity": ["sin", 0, 0.1, 5, 0, 0, 0, 0, 0, 1],
             "velocity": ["file", "./study_fish/bodyC.dat"],
             "objfile": objfolder + "/bodyC20.obj"}

    inputfiles = [waterA, waterB, waterC, bodyA, bodyB, bodyC]
    # inputfiles = [waterA, waterC, bodyA, bodyB, bodyC]

    setting = {"max_dt": 0.01,
               "end_time_step": 10000,
               "end_time": 9}

    generate_input_files(inputfiles, setting, IO_dir, id)
elif "fish" in SimulationCase:
    start = 0.

    id = SimulationCase

    objfolder = code_home_dir + "/cpp/obj/fish"
    water = {"name": "water",
             "type": "Fluid",
             "objfile": objfolder + "/water100.obj"}

    tank = {"name": "tank",
            "type": "RigidBody", "isFixed": True,
            "objfile": objfolder + "/tank50.obj"}

    bodyA = {"name": "bodyA",
             "type": "RigidBody",
             "COM": [0., 0, 0.25],
             "mass": 10**10,
             "MOI": [10**10, 10**10, 10**10],
             "output": "json",
             #  "velocity": ["sin", 0, 0.1, 5, 0, 0, 0, 0, 0, 1],
             "velocity": ["file", "./study_fish/bodyA.dat"],
             "objfile": objfolder + "/bodyA50.obj"}

    bodyB = {"name": "bodyB",
             "type": "RigidBody",
             "COM": [0.35, 0, 0.25],
             "mass": 10**10,
             "MOI": [10**10, 10**10, 10**10],
             "output": "json",
             #  "velocity": ["sin", 0, 0.1, 5, 0, 0, 0, 0, 0, 1],
             "velocity": ["file", "./study_fish/bodyB.dat"],
             "objfile": objfolder + "/bodyB50.obj"}

    bodyC = {"name": "bodyC",
             "type": "RigidBody",
             "COM": [0.7, 0, 0.25],
             "mass": 10**10,
             "MOI": [10**10, 10**10, 10**10],
             "output": "json",
             #  "velocity": ["sin", 0, 0.1, 5, 0, 0, 0, 0, 0, 1],
             "velocity": ["file", "./study_fish/bodyC.dat"],
             "objfile": objfolder + "/bodyC50.obj"}

    inputfiles = [tank, water, bodyA, bodyB, bodyC]

    setting = {"max_dt": 0.01,
               "end_time_step": 10000,
               "end_time": 9}

    generate_input_files(inputfiles, setting, IO_dir, id)
elif "testALE" in SimulationCase:
    objfolder = code_home_dir + "/cpp/obj/testALE"

    water = {"name": "water",
             "type": "Fluid",
             "objfile": objfolder + "/water500mod.obj"}

    tank = {"name": "tank",
            "type": "RigidBody",
            "velocity": ["const", 0, 0.0, 0, 0, 0, 0, 0, 1],
            "objfile": objfolder + "/tank100.obj"}

    cylinder = {"name": "cylinder",
                "type": "RigidBody",
                # "velocity": ["sin", 0, 0.05, 5, 0, 1, 0],
                "velocity": ["const", 0, pi/180., 0, 0, 0, 0, 0, 1],
                "objfile": objfolder + "/cylinder100.obj"}

    cuboid = {"name": "cuboid",
              "type": "RigidBody",
              #   "velocity": ["sin", 0, -0.05, 5, 0, 1, 0],
              "velocity": ["const", 0, pi/180., 0, 0, 0, 0, 0, 1],
              "objfile": objfolder + "/cuboid100.obj"}

    id = SimulationCase

    inputfiles = [water, tank, cylinder, cuboid]

    setting = {"max_dt": 0.02,
               "end_time_step": 10000,
               "end_time": 9}

    generate_input_files(inputfiles, setting, IO_dir, id)
elif "simple_barge" in SimulationCase:
    input_dir += SimulationCase
    os.makedirs(input_dir, exist_ok=True)
    output_dir = sys_home_dir + "/BEM/" + SimulationCase
    os.makedirs(output_dir, exist_ok=True)
    id = SimulationCase
    id = add_id(id)
    # id += '_with_mooring'
    water = {"name": "water", "type": "Fluid"}
    tank = {"name": "tank", "type": "RigidBody", "isFixed": True}

    start = 0.
    a = 1.
    T = 7.
    h = 80
    z_surface = 80
    wavemaker = {"name": "wavemaker",
                 "type": "RigidBody",
                 "velocity": ["piston", start, a, T, h, 1, 0, 0]
                 }

    # 係留索の設定

    X_anchorA = [50., 50., 0.]
    X_fair_leadA = [75., 50., 71.]

    X_anchorB = [150., 50., 0.]
    X_fair_leadB = [125., 50, 71.]

    stiffness = 14 * 10**8  # ! [N/m]
    damp = .5  # ! [N/(m/s^2)]
    density = 100.  # ! [kg/m]
    total_lengthA = math.sqrt((X_anchorA[0]-X_fair_leadA[0])**2 + (
        X_anchorA[1]-X_fair_leadA[1])**2 + (X_anchorA[2]-X_fair_leadA[2])**2)
    total_lengthB = math.sqrt((X_anchorB[0]-X_fair_leadB[0])**2 + (
        X_anchorB[1]-X_fair_leadB[1])**2 + (X_anchorB[2]-X_fair_leadB[2])**2)
    n_points = 30
    diam = 0.1

    if "_with_mooring" in id:
        floatingbody = {"name": "float",
                        "type": "RigidBody",
                        "velocity": "floating",
                        "mooringA": ["mooringA",
                                     X_anchorA[0], X_anchorA[1], X_anchorA[2],
                                     X_fair_leadA[0], X_fair_leadA[1], X_fair_leadA[2],
                                     total_lengthA,
                                     n_points,
                                     density,
                                     stiffness,
                                     damp,
                                     diam],
                        "mooringB": ["mooringB",
                                     X_anchorB[0], X_anchorB[1], X_anchorB[2],
                                     X_fair_leadB[0], X_fair_leadB[1], X_fair_leadB[2],
                                     total_lengthB,
                                     n_points,
                                     density,
                                     stiffness,
                                     damp,
                                     diam],
                        }
    else:
        floatingbody = {"name": "float",
                        "type": "RigidBody",
                        "velocity": "floating"}

    A = 2450.00
    floatingbody["mass"] = m = (rho*g*7.5*A)/g
    floatingbody["COM"] = [100., 50., 75.]
    floatingbody["radius_of_gyration"] = [20., 20., 20.]
    floatingbody["MOI"] = [m*math.pow(floatingbody["radius_of_gyration"][0], 2),
                           m*math.pow(floatingbody["radius_of_gyration"][1], 2),
                           m*math.pow(floatingbody["radius_of_gyration"][2], 2)]

    objfolder = code_home_dir + "/cpp/obj/tsukada2022_no_pool_small_case"
    water["objfile"] = objfolder + "/water100_modmod.obj"
    wavemaker["objfile"] = objfolder + "/wavemaker20.obj"
    tank["objfile"] = objfolder + "/wavetank10.obj"
    floatingbody["objfile"] = objfolder+"/floatingbody10.obj"

    inputfiles = [tank, wavemaker, water, floatingbody]

    setting = {"max_dt": dt,
               "end_time_step": 100000,
               "end_time": 100,
               "ALEPERIOD": ALEPERIOD}

    generate_input_files(inputfiles, setting, IO_dir, id)
elif "moon_pool" in SimulationCase:
    for T in [5 + 0.5 * i for i in range(0, 11)]:

        input_dir = "./input_files/"
        start = 0.
        a = 0.8
        h = 80
        z_surface = 80

        # pool_size = "large"
        pool_size = "none"

        if pool_size == "large":
            id = "_large"
        elif pool_size == "none":
            id = "_no"

        id += "_a" + str(a).replace(".", "d")
        id += "_T" + str(T).replace(".", "d")
        id += "_h" + str(h)

        # id += "_modified_mesh"

        input_dir += SimulationCase + id
        os.makedirs(input_dir, exist_ok=True)
        output_dir = sys_home_dir + "/BEM/" + SimulationCase + id
        os.makedirs(output_dir, exist_ok=True)

        water = {"name": "water", "type": "Fluid"}

        tank = {"name": "tank", "type": "RigidBody", "isFixed": True}

        wavemaker = {"name": "wavemaker",
                     "type": "SoftBody",
                     "isFixed": True,
                     # "isFixed": True}
                     "velocity": ["linear_traveling_wave", start, a, T, h, z_surface]}

        floatingbody = {"name": "float",
                        "type": "RigidBody",
                        # "velocity": ["sin", 0, a, T]}
                        "velocity": "floating"}

        # 浮体の種類
        if pool_size == "large":
            objfolder = code_home_dir + "/cpp/obj/tsukada2022_large_pool"
            water["objfile"] = objfolder + "/water1000mod.obj"
            wavemaker["objfile"] = objfolder + "/wavemaker100.obj"
            tank["objfile"] = objfolder + "/tank10.obj"
            floatingbody["objfile"] = objfolder+"/floating_body50.obj"
            A = 1528.00
            floatingbody["mass"] = m = (1000.*g*7.5*A)/g
            floatingbody["COM"] = [200., 75., 75.]
            floatingbody["radius_of_gyration"] = [20., 20., 20.]
        elif pool_size == "none" or pool_size == "no":
            objfolder = code_home_dir + "/cpp/obj/tsukada2022_no_pool"
            water["objfile"] = objfolder + "/water960mod.obj"
            wavemaker["objfile"] = objfolder + "/wavemaker100.obj"
            tank["objfile"] = objfolder + "/tank10.obj"
            floatingbody["objfile"] = objfolder+"/floating_body_310.obj"
            A = 2450.00
            floatingbody["mass"] = m = (rho*g*7.5*A)/g
            floatingbody["COM"] = [200., 75., 75.]
            floatingbody["radius_of_gyration"] = [20., 20., 20.]

        floatingbody["MOI"] = [m*math.pow(floatingbody["radius_of_gyration"][0], 2),
                               m*math.pow(floatingbody["radius_of_gyration"][1], 2),
                               m*math.pow(floatingbody["radius_of_gyration"][2], 2)]

        inputfiles = [tank, wavemaker, water, floatingbody]

        setting = {"max_dt": dt,
                   "end_time_step": 1000,
                   "end_time": 100,
                   "ALEPERIOD": 1}

        id = SimulationCase
        generate_input_files(inputfiles, setting, IO_dir, id)
elif "two_floatingbodies" in SimulationCase:
    start = 0.
    a = 1.5
    T = 7.  # 5-8
    h = 150
    z_surface = 150

    id = "_a" + str(a).replace(".", "d")\
        + "_T" + str(T).replace(".", "d")\
        + "_h" + str(h).replace(".", "d")

    input_dir += SimulationCase + id
    os.makedirs(input_dir, exist_ok=True)
    output_dir = sys_home_dir + "/BEM/" + SimulationCase + id
    os.makedirs(output_dir, exist_ok=True)

    water = {"name": "water",
             "type": "Fluid"}

    tank = {"name": "tank",
            "type": "RigidBody",
            "isFixed": True}

    wavemaker = {"name": "wavemaker",
                 "type": "SoftBody",
                 "isFixed": True,
                 "velocity": ["linear_traveling_wave", start, a, T, h, z_surface]}

    floatingbody_a = {"name": "float_a",
                      "COM": [300., 150., 150-78./2]}

    floatingbody_b = {"name": "float_b",
                      "COM": [600., 150., 150-78./2]}

    floatingbody_c = {"name": "float_c",
                      "COM": [900., 150., 150-78./2]}

    A = 170.779
    for x in [floatingbody_a, floatingbody_b, floatingbody_c]:
        x["type"] = "RigidBody"
        x["velocity"] = "floating"
        x["mass"] = m = (1000.*g*78*A)/g
        x["radius_of_gyration"] = [20., 20., 20.]
        x["MOI"] = [m*math.pow(x["radius_of_gyration"][0], 2),
                    m*math.pow(x["radius_of_gyration"][1], 2),
                    m*math.pow(x["radius_of_gyration"][2], 2)]

    # ------------------------------------------------------------------------#

    objfolder = code_home_dir + "/cpp/obj/2022Tonegawa/two_floatingbodies/"
    water["objfile"] = objfolder + "/water300_two_mod2.obj"
    wavemaker["objfile"] = objfolder + "/wavemaker100_two.obj"
    tank["objfile"] = objfolder + "/tank10_two.obj"
    floatingbody_a["objfile"] = objfolder + "floatingbody_a50.obj"
    floatingbody_b["objfile"] = objfolder + "floatingbody_b50_two.obj"
    floatingbody_c["objfile"] = objfolder + "floatingbody_c50_two.obj"

    inputfiles = [tank, wavemaker, water,
                  floatingbody_b, floatingbody_c]

    setting = {"max_dt": 0.3,
               "end_time_step": 12000,
               "end_time": 120}

    id = SimulationCase
    generate_input_files(inputfiles, setting, IO_dir, id)
elif "Gu2018Float1_offset_neg0d074" in SimulationCase:
    input_dir += SimulationCase
    os.makedirs(input_dir, exist_ok=True)
    output_dir = sys_home_dir + "/BEM/" + SimulationCase
    os.makedirs(output_dir, exist_ok=True)

    water = {"name": "water", "type": "Fluid"}
    tank = {"name": "tank", "type": "RigidBody", "isFixed": True}

    float = {"name": "float",
             "type": "RigidBody",
             "velocity": "floating"}

    float["mass"] = m = 9.75
    # float["COM"] = [0., 0., 0.09+0.8]  # 今回は重要ではない
    float["COM"] = [0., 0., 0.]  # 今回は重要ではない
    float["radius_of_gyration"] = [10**10, 10**10, 10**10]
    float["MOI"] = [m*math.pow(float["radius_of_gyration"][0], 2),
                    m*math.pow(float["radius_of_gyration"][1], 2),
                    m*math.pow(float["radius_of_gyration"][2], 2)]

    objfolder = code_home_dir + "/cpp/obj/" + SimulationCase
    water["objfile"] = objfolder + "/water300.obj"
    tank["objfile"] = objfolder + "/tank10.obj"
    float["objfile"] = objfolder+"/float10.obj"

    inputfiles = [tank, water, float]

    setting = {"max_dt": 0.02,
               "end_time_step": 10000,
               "end_time": 4}
    id = SimulationCase
    generate_input_files(inputfiles, setting, IO_dir, id)
elif "Tonegawa2022" in SimulationCase:
    input_dir += SimulationCase
    os.makedirs(input_dir, exist_ok=True)
    output_dir = sys_home_dir + "/BEM/2022Tonegawa3"
    os.makedirs(output_dir, exist_ok=True)

    water = {"name": "water",
             "type": "Fluid"}

    tank = {"name": "tank",
            "type": "RigidBody",
            "isFixed": True}

    start = 0.
    a = 1.
    T = 7.
    h = 150
    z_surface = 150

    wavemaker = {"name": "wavemaker",
                 "type": "SoftBody",
                 "isFixed": True,
                 "velocity": ["linear_traveling_wave", start, a, T, h, z_surface]}

    floatingbody = {"name": "float",
                    "type": "RigidBody",
                    "velocity": "floating"}

    A = 162.122
    floatingbody["mass"] = m = (1000.*g*78*A)/g
    floatingbody["COM"] = [75., 50., 150-78./2]
    floatingbody["radius_of_gyration"] = [20., 20., 20.]
    floatingbody["MOI"] = [m*math.pow(floatingbody["radius_of_gyration"][0], 2),
                           m*math.pow(floatingbody["radius_of_gyration"][1], 2),
                           m*math.pow(floatingbody["radius_of_gyration"][2], 2)]

    objfolder = code_home_dir + "/cpp/obj/2022Tonegawa/three_d100"
    water["objfile"] = objfolder + "/water150_mod.obj"
    wavemaker["objfile"] = objfolder + "/wavemaker100.obj"
    tank["objfile"] = objfolder + "/tank10.obj"
    floatingbody["objfile"] = objfolder+"/floatingbody50.obj"

    inputfiles = [tank, wavemaker, water, floatingbody]

    setting = {"max_dt": 0.2,
               "end_time_step": 10000,
               "end_time": 10000}
    id = SimulationCase
    generate_input_files(inputfiles, setting, IO_dir, id)
elif "three_floatingbodies" in SimulationCase:
    start = 0.
    a = 1.5
    T = 6.0  # 5-8
    h = 150
    z_surface = 150

    id = "_a" + str(a).replace(".", "d")\
        + "_T" + str(T).replace(".", "d")\
        + "_h" + str(h).replace(".", "d")

    input_dir += SimulationCase + id
    os.makedirs(input_dir, exist_ok=True)
    output_dir = sys_home_dir + "/BEM/" + SimulationCase + id
    os.makedirs(output_dir, exist_ok=True)

    water = {"name": "water",
             "type": "Fluid"}

    tank = {"name": "tank",
            "type": "RigidBody",
            "isFixed": True}

    wavemaker = {"name": "wavemaker",
                 "type": "SoftBody",
                 "isFixed": True,
                 "velocity": ["linear_traveling_wave", start, a, T, h, z_surface]}

    floatingbody_a = {"name": "float_a",
                      "COM": [300., 150., 150-78./2]}

    floatingbody_b = {"name": "float_b",
                      "COM": [600., 150., 150-78./2]}

    floatingbody_c = {"name": "float_c",
                      "COM": [900., 150., 150-78./2]}

    A = 170.779
    for x in [floatingbody_a, floatingbody_b, floatingbody_c]:
        x["type"] = "RigidBody"
        x["velocity"] = "floating"
        x["mass"] = m = (1000.*g*78*A)/g
        x["radius_of_gyration"] = [20., 20., 20.]
        x["MOI"] = [m*math.pow(x["radius_of_gyration"][0], 2),
                    m*math.pow(x["radius_of_gyration"][1], 2),
                    m*math.pow(x["radius_of_gyration"][2], 2)]

    # ------------------------------------------------------------------------#

    objfolder = code_home_dir + "/cpp/obj/2022Tonegawa/three_d100/"
    water["objfile"] = objfolder + "/water0_200.obj"
    wavemaker["objfile"] = objfolder + "/wavemaker100.obj"
    tank["objfile"] = objfolder + "/tank10.obj"
    floatingbody_a["objfile"] = objfolder + "floatingbody_a50.obj"
    floatingbody_b["objfile"] = objfolder + "floatingbody_b50.obj"
    floatingbody_c["objfile"] = objfolder + "floatingbody_c50.obj"

    inputfiles = [tank, wavemaker, water,
                  floatingbody_a, floatingbody_b, floatingbody_c]

    setting = {"max_dt": 0.3,
               "end_time_step": 12000,
               "end_time": 150}

    id = SimulationCase
    generate_input_files(inputfiles, setting, IO_dir, id)
elif "Retzler2000simple" in SimulationCase:
    input_dir += SimulationCase
    os.makedirs(input_dir, exist_ok=True)
    output_dir = sys_home_dir + "/BEM/" + SimulationCase
    os.makedirs(output_dir, exist_ok=True)

    water = {"name": "water",
             "type": "Fluid",
             }

    tank = {"name": "tank",
            "type": "RigidBody",
            "isFixed": True}

    wavemaker = {"name": "wavemaker",
                 "type": "RigidBody",
                 #  "isFixed": True,
                 "velocity": ["Retzler2000", 0.1]}

    objfolder = code_home_dir + "/cpp/obj/chaplin2000simple/"
    water["objfile"] = objfolder+"/water_refined.obj"
    wavemaker["objfile"] = objfolder+"/cylinder200.obj"
    tank["objfile"] = objfolder+"/tank.obj"
    wavetank_ignore = False
    inputfiles = [tank, wavemaker, water]
    setting = {"max_dt": dt,
               "end_time_step": 10000,
               "end_time": 0.5,
               "ALEPERIOD": ALEPERIOD}

    id = SimulationCase
    id = add_id(id)
    generate_input_files(inputfiles, setting, IO_dir, id)
elif "Retzler2000" in SimulationCase:
    input_dir += SimulationCase
    os.makedirs(input_dir, exist_ok=True)
    output_dir = sys_home_dir + "/BEM/Retzler2000simple"
    os.makedirs(output_dir, exist_ok=True)

    water = {"name": "water",
             "type": "Fluid",
             }

    tank = {"name": "tank",
            "type": "RigidBody",
            "isFixed": True}

    wavemaker = {"name": "wavemaker",
                 "type": "RigidBody",
                 #  "isFixed": True,
                 "velocity": ["Retzler2000", 0.02]}

    objfolder = code_home_dir + "/cpp/obj/chaplin2000/"
    water["objfile"] = objfolder+"/water_remeshed2.obj"
    wavemaker["objfile"] = objfolder+"/cylinder200.obj"
    tank["objfile"] = objfolder+"/tank.obj"
    wavetank_ignore = False
    inputfiles = [tank, wavemaker, water]
    setting = {"max_dt": 0.002,
               "end_time_step": 100,
               "end_time": 0.4}
    id = SimulationCase
    generate_input_files(inputfiles, setting, IO_dir, id)
elif "three_200" in SimulationCase:
    start = 0.
    a = 1.5
    T = 6.0  # 5-8
    h = 150
    z_surface = 150

    id = "_a" + str(a).replace(".", "d")\
        + "_T" + str(T).replace(".", "d")\
        + "_h" + str(h).replace(".", "d")\
        + "_"

    input_dir += SimulationCase + id
    os.makedirs(input_dir, exist_ok=True)
    output_dir = sys_home_dir + "/BEM/" + SimulationCase + id
    os.makedirs(output_dir, exist_ok=True)

    water = {"name": "water",
             "type": "Fluid"}

    tank = {"name": "tank",
            "type": "RigidBody",
            "isFixed": True}

    wavemaker = {"name": "wavemaker",
                 "type": "SoftBody",
                 "isFixed": True,
                 "velocity": ["linear_traveling_wave", start, a, T, h, z_surface]}

    floatingbody_a = {"name": "floatingbody_a",
                      "COM": [400., 150., 150-78./2]}

    floatingbody_b = {"name": "floatingbody_b",
                      "COM": [600., 150., 150-78./2]}

    floatingbody_c = {"name": "floatingbody_c",
                      "COM": [800., 150., 150-78./2]}

    A = 170.779
    for x in [floatingbody_a, floatingbody_b, floatingbody_c]:
        x["type"] = "RigidBody"
        x["velocity"] = "floating"
        x["mass"] = m = (1000.*g*78*A)/g
        x["radius_of_gyration"] = [20., 20., 20.]
        x["MOI"] = [m*math.pow(x["radius_of_gyration"][0], 2),
                    m*math.pow(x["radius_of_gyration"][1], 2),
                    m*math.pow(x["radius_of_gyration"][2], 2)]

    # ------------------------------------------------------------------------#

    objfolder = code_home_dir + "/cpp/obj/2022Tonegawa/three_d200/"
    water["objfile"] = objfolder + "/water_modified2.obj"
    wavemaker["objfile"] = objfolder + "/wavemaker100.obj"
    tank["objfile"] = objfolder + "/tank10_200.obj"
    floatingbody_a["objfile"] = objfolder + "floatingbodyA50_200.obj"
    floatingbody_b["objfile"] = objfolder + "floatingbodyB50_200.obj"
    floatingbody_c["objfile"] = objfolder + "floatingbodyC50_200.obj"

    inputfiles = [tank, wavemaker, water,
                  floatingbody_a, floatingbody_b, floatingbody_c]

    setting = {"max_dt": 0.5,
               "end_time_step": 12000,
               "end_time": 120}

    id = SimulationCase
    generate_input_files(inputfiles, setting, IO_dir, id)
elif "Goring1979" in SimulationCase:

    r'''DOC_EXTRACT 2_2_0_validation_Goring1979

    ## Goring1979

    The case is based on Goring's thsis \cite{Goring1979}.
    We simulated the case to see the efficiency of the pseudo-quad elements against the linear elements.
    Different element types can be used for BIE and ALE separately.
    
    The cases are as follows:

    | the upstream length | the downstream length | the upstream depth $`h_1`$ | the height of the shelf | Gauge1 | Gauge2 | Gauge3 |
    |:---:|:---:|:---:|:---:|:---:|:---:|:---:|
    | 13 m | 19.84 m | 25 cm | 15.54 cm | 23 $`h_1`$ (5.75 m) | at the step | 60 $`h_1`$ (5.68 m) |

    <p align="center">
    <img src="./img/Goring1979Fig5.1p122.png" height="400">
    <img src="./img/Goring1979Fig5.2p123.png" height="400">
    </p>

    The downstream length of the shelf is too long, 19.84 m, therefore, in this simulation case, the length is reduced to 10 m.

    [CAD data](https://a360.co/3EZSCBN)
    
    inputファイルの生成

    ```sh
    python3.11 input_generator.py Goring1979 -m water0d1.obj -dt 0.05 -e pseudo_quad -o /Volumes/home/BEM/Goring1979/
    python3.11 input_generator.py Goring1979 -m water0d1.obj -dt 0.05 -o /Volumes/home/BEM/Goring1979/
    python3.11 input_generator.py Goring1979 -m water0d09.obj -dt 0.05 -o /Volumes/home/BEM/Goring1979/
    python3.11 input_generator.py Goring1979 -m water0d09.obj -dt 0.05 -o /Volumes/home/BEM/Goring1979/
    ```

    実行ファイルが`fast`の場合，以下のように実行する．

    ```sh
    ./fast ./input_files/Goring1979_DT0d05_MESHwater0d09.obj_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
    ```    

    '''

    # IDの生成
    id = SimulationCase
    id = add_id(id)
    max_dt = dt

    objfolder = code_home_dir + "/cpp/obj/Goring1979/"

    water = {"name": "water",
             "type": "Fluid",
             "objfile": objfolder + "/water.obj"}

    if meshname != "":
        water["objfile"] = objfolder + "/" + meshname

    vertices, triangles = read_obj(water["objfile"])
    water["vertices"] = vertices
    water["triangles"] = triangles

    tank = {"name": "tank",
            "type": "RigidBody",
            "isFixed": True,
            "objfile": objfolder + "/tank.obj"}


    h = 0.25
    wavemaker = {"name": "wavemaker",
                 "type": "RigidBody",
                 "velocity": ["Goring1979", 3., 0.1*h, h],
                 "objfile": f"{objfolder}/wavemaker.obj"}

    gauges = []
    gauges.append({"name": "gauge1", 
                   "type": "wave gauge",
                  "position": [-5.75, 0, 0.4, -5.75, 0, 0.2]})
    gauges.append({"name": "gauge2",                    
                   "type": "wave gauge",
                  "position": [0, 0, 0.4, 0, 0, 0.2]})
    gauges.append({"name": "gauge3",
                   "type": "wave gauge",
                   "position": [5.68, 0, 0.4, 5.68, 0, 0.2]})

    inputfiles = [tank, wavemaker, water] + gauges

    setting = {"max_dt": max_dt,
               "end_time_step": 100000,
               "end_time": 20,
               "element": element,
               "ALE": ALE,
               "ALEPERIOD": ALEPERIOD}

    generate_input_files(inputfiles, setting, IO_dir, id)
elif "Horikawa2024" in SimulationCase:

    r'''DOC_EXTRACT 2_1_0_input_generator

    ## Horikawa2024

    ### inputファイルの生成

    ```sh
    python3 input_generator.py Horikawa2024 -dt 0.05 -ALEPERIOD 1 -ALE linear -element linear -output ~/BEM/Horikawa2024
    ```

    または，省略して次のように実行する．

    ```sh
    python3 input_generator.py Horikawa2024 -dt 0.05 -output ~/BEM/Horikawa2024
    ```

    <img src="how_to_run_example.png" width="400">

    `/input_files/Horikawa2024_a0d003_T0d625_DT0d05_ELEMlinear_ALElinear_ALEPERIOD1`をコピーして，実行時に指定する．

    ### 実行ファイルの生成（毎回する必要はない）

    `sh clean`で古いファイルを削除する．`cmake`を使ってcppのコンパイル方法を設定する．
    cmakeをした後に，コンパイル（make）する．以下の場合，実行ファイル名はhorikawaとなる．

    ```sh
    sh clean
    cmake -DCMAKE_BUILD_TYPE=Release ../ -DOUTPUT_NAME=horikawa
    make
    ```

    実行ファイル`horikawa`を使って事項する．ただし，`input_files`ディレクトリに保存されているinputファイルを指定して実行する．

    ### 実行

    ```sh
    ./horikawa ./input_files/Horikawa2024_a0d003_T0d625_DT0d05_ELEMlinear_ALElinear_ALEPERIOD1
    ```

    '''

    objfolder = code_home_dir + "/cpp/obj/Horikawa2024"

    # 水とタンク、浮体の設定
    water = {"name": "water",
             "type": "Fluid",
             "objfile": objfolder + "/water.obj"}
    # "objfile": objfolder + "/waterHigh.obj"}

    tank = {"name": "tank",
            "type": "RigidBody",
            "isFixed": True,
            "objfile": objfolder + "/tank.obj"}

    start = 0.
    a = 0.003
    T = 0.625
    float = {"name": "float",
             "type": "RigidBody",
             "velocity": ["sin", start, a, T, 1, 0, 0, 0, 0, 0],
             "objfile": objfolder + "/float.obj"}

    # IDの生成
    id = SimulationCase
    id += "_a" + str(a).replace(".", "d")
    id += "_T" + str(T).replace(".", "d")
    id = add_id(id)

    # 入力ファイルの設定
    inputfiles = [tank, water, float]

    # シミュレーション設定
    setting = {"max_dt": dt,  # デフォルトの時間刻み幅
               "end_time_step": 10000,
               "end_time": 10,
               "element": element,
               "ALE": ALE,
               "ALEPERIOD": ALEPERIOD}

    # 入力ファイルを生成
    generate_input_files(inputfiles, setting, IO_dir, id)
elif "Tonegawa2024Experiment" in SimulationCase:

    r'''DOC_EXTRACT 2_1_0_input_generator

    ## Tonegawa2024

    ### 浮体

    | 項目               | 値                          |
    |--------------------|-----------------------------|
    | サイズ (WxLxH)     | 0.375 x 0.375 x 0.84 m      |
    | ムーンプールサイズ (WxL) | 0.25 x 0.25 m            |
    | 喫水               | 0.056 m                     |
    | 重量               | 7.56 kg (0d075) / 7.0 kg (0d125) / 4.375 kg (0d25) |
    | 慣性モーメント     | [0.1139241482, 0.1139241482, 0.2117076729] (0d075) / [0.1078722037, 0.1078722037, 0.2013731481] (0d125) / [0.07943494545, 0.07943494545, 0.1514044693] (0d25) |
    | 重心位置           | (x, y, z) = (2, 0, h - draft + H_float / 2) |

    ### 水槽

    | 項目               | 値                          |
    |--------------------|-----------------------------|
    | サイズ (WxLxDepth) | 2.15 x 4.45 x 0.6 m         |

    ### 使用するobjファイル

    | ファイル名                       |
    |----------------------------------|
    | Tonegawa2024/water0d075.obj (0d075) / Tonegawa2024/water0d125.obj (0d125) / Tonegawa2024/water0d25.obj (0d25) |
    | Tonegawa2024/float0d075.obj (0d075) / Tonegawa2024/float0d125.obj (0d125) / Tonegawa2024/float0d25.obj (0d25) |
    | Tonegawa2024/tank.obj            |
    | Tonegawa2024/wavemaker.obj       |

    ### 造波装置の運動

    - (x, z) = (-0.31, 0)を通るy軸を中心として，wavemakerの運動をflap型造波で行う．

    '''

    start = 0.
    H = 0.015  # 波高
    # if not degined
    try:
        T
    except:
        T = 0.6   # 波周期

    a = H / 2  # 振幅
    h = 0.6   # 水深
    z_surface = h

    # 浮体設定
    W, L, H_float = 0.375, 0.375, 0.084
    draft = 0.056
    # mass = 4.375
    # mass = 7.
    # COM = [2, 0, h - draft + H_float / 2]
    COM = [1.8675, 0, h - draft + H_float / 2]
    # MOI = [0.07943494545, 0.07943494545, 0.1514044693]
    # MOI = [0.1078722037,0.1078722037,0.2013731481]

    # 水槽設定
    tank_dim = [2.15, 4.45, 0.6]

    # obj ファイル
    objfolder = code_home_dir + "/cpp/obj/Tonegawa2024"

    if "0d075" in suffix:
        # 面積	5642.64 cm^2
        # 密度	1.493 g / cm^3
        # 質量	7559.476 g
        # 体積	5064.84 cm^3
        # 重心	0.00 cm, 0.00 cm, 2.577 cm
        # 重心の慣性モーメント   (g cm^2)
        # 	Ixx = 1.074E+06
        # 	Ixy = -6.950E-10
        # 	Ixz = 2.289E-10
        # 	Iyx = -6.950E-10
        # 	Iyy = 1.074E+06
        # 	Iyz = -4.778E-10
        # 	Izx = 2.289E-10
        # 	Izy = -4.778E-10
        # 	Izz = 2.063E+06
        float_obj = objfolder + "/float0d075.obj"
        water_obj = objfolder + "/water0d075.obj"
        mass = 7.56
        # MOI = [0.0899209,0.0899209,0.163718]
        MOI = [1.074 * 1e6, 1.074 * 1e6, 2.063 * 1e6]
        MOI[0] *= 1e-7
        MOI[1] *= 1e-7
        MOI[2] *= 1e-7
        COM[2] = h - draft + 2.577 / 100
    if "0d125" in suffix:
        # 質量	7002.998 g
        # 体積	4692.00 cm^3
        # 重心	0.00 cm, 0.00 cm, 2.744 cm
        # 重心の慣性モーメント   (g cm^2)
        #     Ixx = 1.040E+06
        #     Ixy = -1.390E-09
        #     Ixz = 1.910E-10
        #     Iyx = -1.390E-09
        #     Iyy = 1.040E+06
        #     Iyz = -3.909E-10
        #     Izx = 1.910E-10
        #     Izy = -3.909E-10
        #     Izz = 1.991E+06
        float_obj = objfolder + "/float0d125.obj"
        water_obj = objfolder + "/water0d125.obj"
        mass = 7.
        # MOI = [0.0899375,0.0899375,0.16452]
        MOI = [1.040 * 1e6, 1.040 * 1e6, 1.991 * 1e6]
        MOI[0] *= 1e-7
        MOI[1] *= 1e-7
        MOI[2] *= 1e-7
        COM[2] = h - draft + 2.744 / 100
    elif "0d25" in suffix:
        # 質量 4375.008 g
        # 体積 2931.25 cm^3
        # 重心 0.00 cm, 0.00 cm, 4.021 cm
        # 重心の慣性モーメント (g cm^2)
        # Ixx = 8.040 E + 05
        # Ixy = -1.390 E - 09
        # Ixz = 0.00
        # Iyx = -1.390 E - 09
        # Iyy = 8.040 E + 05
        # Iyz = -3.909 E - 10
        # Izx = 0.00
        # Izy = -3.909 E - 10
        # Izz = 1.529 E + 06
        float_obj = objfolder + "/float0d25.obj"
        water_obj = objfolder + "/water0d25.obj"
        mass = 4.375
        # MOI = [0.0985854,0.0985854,0.182528]
        MOI = [8.040 * 1e5, 8.040 * 1e5, 1.529 * 1e6]
        MOI[0] *= 1e-7
        MOI[1] *= 1e-7
        MOI[2] *= 1e-7
        COM[2] = h - draft + 4.021 / 100

    water = {"name": "water",
             "type": "Fluid",
             "objfile": water_obj}

    float = {"name": "float",
             "type": "RigidBody",
             "velocity": "floating",
             "mass": mass,
             "COM": COM,
             "MOI": MOI,
             #  "isFixed": [True,True,False,False,False,False],
             "objfile": float_obj}

    if "linear_cable" in suffix:
        Wf = 0.375
        k = 10*0.5/85*9.81
        x = float["COM"][0]
        y = float["COM"][1]
        z = float["COM"][2]
        float["linear_cable1"] = [x + Wf/2, y +
                                  Wf/2, z, x + Wf/2 + 2, y + Wf/2, 0, k]
        float["linear_cable2"] = [x - Wf/2, y +
                                  Wf/2, z, x - Wf/2 - 2, y + Wf/2, 0, k]
        float["linear_cable3"] = [x + Wf/2, y -
                                  Wf/2, z, x + Wf/2 + 2, y - Wf/2, 0, k]
        float["linear_cable4"] = [x - Wf/2, y -
                                  Wf/2, z, x - Wf/2 - 2, y - Wf/2, 0, k]

    tank = {"name": "tank", "type": "RigidBody",
            "isFixed": True, "objfile": objfolder + "/tank.obj"}

    wavemaker = {"name": "wavemaker",
                 "type": "RigidBody",
                 "velocity": ["flap", start, a, T, h, h, 0, 1, 0],
                 "objfile": objfolder + "/wavemaker.obj"}

    bottom_in_z = 0
    h = 0.6
    absorber = {"name": "absorber",
                "type": "Absorber",
                "isFixed": True,
                "wave_theory": [0, T, h, bottom_in_z],
                "objfile": objfolder+"/absorber.obj"}

    gauges = []
    for i in range(8):
        gauges.append({"name": f"gauge{i}",
                       "type": "wave gauge",
                       "position": [0.25*(i+1), 0., 1., 0.25*(i+1), 0., 0.2]})

    inputfiles = [tank, wavemaker, water, float, absorber]
    inputfiles += gauges

    id = SimulationCase
    id += "_T" + "{:.2f}".format(T).replace(".", "d")
    id = add_id(id)

    setting = {"max_dt": dt,
               "end_time_step": 1000000,
               "end_time": 12,
               "ALE": ALE,
               "ALEPERIOD": ALEPERIOD}

    generate_input_files(inputfiles, setting, IO_dir, id)
elif "Fredriksen2015" in SimulationCase:

    r'''DOC_EXTRACT 2_1_0_validation_Fredriksen2015

    ## Fredriksen2015

    '''

    # 波の設定
    start = 0.
    H = 0.05
    eps = 1/60.
    # T = 0.9
    h = 1.
    z_surface = h

    # 浮体の設定
    draft = 0.097  # 喫水
    vol = draft * (0.201*0.586)*2.
    M = rho * vol
    COM = [0, 0, h-draft+0.091]  # 重心位置
    ryy = 0.18  # radius of gyration
    Ixx = 500
    Iyy = M*ryy**2
    Izz = 500

    water = {"name": "water",
             "type": "Fluid"}

    tank = {"name": "tank",
            "type": "RigidBody",
            "isFixed": True}

    # wavemaker = {"name": "wavemaker",
    #              "type": "RigidBody",
    #              "velocity": ["flap_wave_steepness", start, eps, T, h, h, 0, 1, 0]}

    z = 0.05 + 1.
    x = 0.201 + 2.
    float = {"name": "float",
             "type": "RigidBody",
             "velocity": "floating",
             "mass": M,
             "COM": COM,
             "MOI": [Ixx, Iyy, Izz],
             "spring1": [-x, 0.25, z, -x - 0.5, 0.25, z, 43.7, 7.5],
             "spring2": [-x, -0.25, z, -x - 0.5, -0.25, z, 42.5, 7.5],
             "spring3": [x, 0., z, x + 0.5, 0., z, 88.2, 15]}

    h = 1
    WS = 1/60  # H/h = 2a/h -> a = eps*h/2
    bottom_in_z = 0.
    absorber1 = {"name": "absorber1",
                 "type": "Absorber",
                 "isFixed": True,
                 "wave_theory": [WS*h/2, T, h, bottom_in_z]}

    absorber2 = {"name": "absorber2",
                 "type": "Absorber",
                 "isFixed": True,
                 "wave_theory": [0, T, h, bottom_in_z]}

    objfolder = code_home_dir + "/cpp/obj/Fredriksen2015"
    water["objfile"] = objfolder + "/water.obj"
    # wavemaker["objfile"] = objfolder + "/wavemaker.obj"
    tank["objfile"] = objfolder + "/tank.obj"
    float["objfile"] = objfolder + "/float.obj"
    absorber1["objfile"] = objfolder + "/absorber1.obj"
    absorber2["objfile"] = objfolder + "/absorber2.obj"

    inputfiles = [tank, water, float, absorber1, absorber2]

    id = SimulationCase
    id += "_WS" + "{:.5f}".format(WS).replace(".", "d")
    id += "_T" + "{:.2f}".format(T).replace(".", "d")
    id = add_id(id)

    setting = {"max_dt": dt,
               "end_time_step": 1000000,
               "end_time": 15,
               "ALE": ALE,
               "element": element,
               "ALEPERIOD": ALEPERIOD}

    generate_input_files(inputfiles, setting, IO_dir, id)
elif "OC6DeepCwind" in SimulationCase:

    start = 0.0

    # Wave parameters (JONSWAP as default)
    H = 7.4  # Significant wave height in meters
    T = 12.0  # Peak wave period in seconds
    gamma = 3.3  # JONSWAP spectrum parameter

    id0 = ""
    id = SimulationCase + id0 + "_H" + str(H).replace(".", "d")
    id += "_T" + str(T).replace(".", "d")

    water = {"name": "water", "type": "Fluid"}

    tank = {"name": "tank",
            "type": "RigidBody",
            "isFixed": True}

    float = {"name": "float",
             "type": "RigidBody",
             "velocity": "floating"}

    # Semisubmersible parameters
    draft = 20.0  # Approximate draft in meters
    float_H = 38.0  # Overall height in meters
    float_W = 60.0  # Width in meters
    float_L = 60.0  # Length in meters

    rho = 1025.0  # Seawater density in kg/m^3
    displacement_volume = 14053.0  # Displacement volume in m^3

    float["mass"] = m = rho * displacement_volume

    kx, ky, kz = 18.7, 20.3, 21.1  # Radius of gyration values (approx.)

    Ixx = m * kx**2
    Iyy = m * ky**2
    Izz = m * kz**2

    float["MOI"] = [Ixx, Iyy, Izz]
    float["COM"] = [0, 0, -7.32]  # Center of mass (Gz value)

    objfolder = code_home_dir + "/cpp/obj/OC6DeepCwind"

    water["objfile"] = objfolder + "/water.obj"
    tank["objfile"] = objfolder + "/tank.obj"
    float["objfile"] = objfolder + "/float.obj"

    # Wave theory (JONSWAP parameters)
    bottom_in_z = -draft
    wave_theory = [H / 2, T, 80.0, bottom_in_z, gamma]

    absorber = {"name": "absorber",
                "type": "Absorber",
                "isFixed": True,
                "objfile": objfolder + "/absorber.obj",
                "wave_theory": wave_theory}

    inputfiles = [tank, water, float, absorber]

    id = add_id(id)

    setting = {"max_dt": dt,
               "end_time_step": 1000000,
               "end_time": 100,
               "element": element,
               "ALE": ALE,
               "ALEPERIOD": ALEPERIOD}

    generate_input_files(inputfiles, setting, IO_dir, id)
elif "DeepCWind" in SimulationCase:

    r'''DOC_EXTRACT_2_1_0_validation_DeepCWind

    ## DeepCWind

    このケースは，\cite{Wang2022}に基づいてる．
    Wang(2022)のTable1によると，

    | 項目 | 値 |
    |:---:|:---:|
    | 全体質量 | 1.4196E+7 kg |
    | 排水量 | 14053 m3 |
    | （水面が原点）重心位置 z | -7.32 m |
    | ピッチ慣性モーメント(Iyy) | 1.2979E+10 kg m2 |
            
    このケースでは，

    | 項目 | 値 |
    |:---:|:---:|
    | 全体質量	1.3958E+7 kg |
    | 喫水 (Draft)	20 m |
    | 排水量	1.3917E+4 m3 |
    | （水面が原点）重心位置 z | -8.07 m |
    | ロール慣性モーメント	1.3947E+10 kg-m2 |
    | ピッチ慣性モーメント	1.5552E+10 kg-m2 |
    | ヨー慣性モーメント	1.3692E+10 kg-m2 |

    水深
    波高 H = 7.4 m
    波長 L = 150 m
    波周期

    ```sh
    python3.11 input_generator.py Goring1979 -m water0d1.obj -dt 0.05 -e pseudo_quad -o /Volumes/home/BEM/Goring1979/
    python3.11 input_generator.py Goring1979 -m water0d1.obj -dt 0.05 -o /Volumes/home/BEM/Goring1979/
    python3.11 input_generator.py Goring1979 -m water0d09.obj -dt 0.05 -o /Volumes/home/BEM/Goring1979/
    python3.11 input_generator.py Goring1979 -m water0d09.obj -dt 0.05 -o /Volumes/home/BEM/Goring1979/
    ```

    実行ファイルが`fast`の場合，以下のように実行する．

    ```sh
    ./fast ./input_files/Goring1979_DT0d05_MESHwater0d09.obj_ELEMlinear_ALEpseudo_quad_ALEPERIOD1
    ```    

    '''

    start = 0.0
    bottom_in_z = 0
    h = 180.0 # 水深
    L = 150.0 # 波長
    theta = 180.0
    wave_theory_L = [H / 2, L, h, bottom_in_z, theta]

    # Mass properties
    displacement_volume = (4297.7 * 3 + 663.661)
    mass = 1000 * displacement_volume  # kg

    # Center of Mass
    COM = [0, 0, h - 8.07]

    # Moments of Inertia
    Ixx = 1.3947e10
    Iyy = 1.5552e10
    Izz = 1.3692e10

    id = SimulationCase

    # Define objects
    water = {"name": "water", "type": "Fluid"}

    tank = {"name": "tank",
            "type": "RigidBody",
            "isFixed": True}

    float = {"name": "float",
             "type": "RigidBody",
             "velocity": "floating",
             "mass": mass,
             "COM": COM,
             "MOI": [Ixx, Iyy, Izz]}

    objfolder = code_home_dir + "/cpp/obj/DeepCWind"

    water["objfile"] = objfolder + "/water.obj"
    tank["objfile"] = objfolder + "/tank.obj"
    float["objfile"] = objfolder + "/float.obj"

    # if suffix contains 9floats, then use water9floats.obj
    # and inclrease the number of floats to 9 (-300,300), (-300,0), (-300,-300), (0,300), (0,0), (0,-300), (300,300), (300,0), (300,-300)
    # (0,0) is the original float position

    absorber = {"name": "absorber",
                "type": "Absorber",
                "isFixed": True,
                "objfile": objfolder + "/absorber.obj",
                "wave_theory_L": wave_theory_L}

    if "9floats" in suffix:

        '''
        y
        |--> x          

        00 <--250-> 01 <--250-> 02   

        03 <--250-> 04 <--250-> 05

        06 <--250-> 07 <--250-> 08
        
        '''

        water["objfile"] = objfolder + "/water9floats.obj"
        
        shift_x = 250
        shift_y = 300

        float00 = copy.copy(float)
        float00["name"] = "float00"
        float00["translation"] = s = [-shift_x, shift_y, 0]
        float00["COM"] = [s[0], s[1], COM[2]]

        float01 = copy.copy(float)
        float01["name"] = "float01"
        float01["translation"] = s = [0, shift_y, 0]
        float01["COM"] = [s[0], s[1], COM[2]]

        float02 = copy.copy(float)
        float02["name"] = "float02"
        float02["translation"] = s = [shift_x, shift_y, 0]
        float02["COM"] = [s[0], s[1], COM[2]]

        float03 = copy.copy(float)
        float03["name"] = "float03"
        float03["translation"] = s = [-shift_x, 0, 0]
        float03["COM"] = [s[0], s[1], COM[2]]

        float04 = copy.copy(float)
        float04["name"] = "float04"
        float04["translation"] = s = [0, 0, 0]
        float04["COM"] = [s[0], s[1], COM[2]]

        float05 = copy.copy(float)
        float05["name"] = "float05"
        float05["translation"] = s = [shift_x, 0, 0]
        float05["COM"] = [s[0], s[1], COM[2]]

        float06 = copy.copy(float)
        float06["name"] = "float06"
        float06["translation"] = s = [-shift_x, -shift_y, 0]
        float06["COM"] = [s[0], s[1], COM[2]]

        float07 = copy.copy(float)
        float07["name"] = "float07"
        float07["translation"] = s = [0, -shift_y, 0]
        float07["COM"] = [s[0], s[1], COM[2]]

        float08 = copy.copy(float)
        float08["name"] = "float08"
        float08["translation"] = s = [shift_x, -shift_y, 0]
        float08["COM"] = [s[0], s[1], COM[2]]

        inputfiles = [tank, water, 
                      float00, float01, float02, 
                      float03,float04, float05, 
                      float06, float07, float08, 
                      absorber]

    if "6floats" in suffix:

        '''
        y
        |--> x          

        00 <--250-> 01

        02 <--250-> 03

        04 <--250-> 05
        
        '''

        water["objfile"] = objfolder + "/water6floats.obj"
        
        shift_x = 250
        shift_y = 300

        float00 = copy.copy(float)
        float00["name"] = "float00"
        float00["translation"] = s = [shift_x, shift_y, 0]
        float00["COM"] = [s[0], s[1], COM[2]]

        float01 = copy.copy(float)
        float01["name"] = "float01"
        float01["translation"] = s = [0, shift_y, 0]
        float01["COM"] = [s[0], s[1], COM[2]]

        float02 = copy.copy(float)
        float02["name"] = "float02"  # 修正点
        float02["translation"] = s = [shift_x, 0, 0]
        float02["COM"] = [s[0], s[1], COM[2]]

        float03 = copy.copy(float)
        float03["name"] = "float03"  # 修正点
        float03["translation"] = s = [0, 0, 0]
        float03["COM"] = [s[0], s[1], COM[2]]

        float04 = copy.copy(float)
        float04["name"] = "float04"  # 修正点
        float04["translation"] = s = [shift_x, -shift_y, 0]
        float04["COM"] = [s[0], s[1], COM[2]]

        float05 = copy.copy(float)
        float05["name"] = "float05"  # 修正点
        float05["translation"] = s = [0, -shift_y, 0]
        float05["COM"] = [s[0], s[1], COM[2]]

        inputfiles = [tank, water, 
                      float00, float01, float02, 
                      float03, float04, float05, 
                      absorber]
        

    elif "4floats" in suffix:

        '''
        y
        |--> x          

        00

        01

        02
        
        03

        '''

        water["objfile"] = objfolder + "/water4floats.obj"

        shift_y = 300

        float00 = copy.copy(float)
        float00["name"] = "float00"
        float00["translation"] = s = [0, shift_y, 0]
        float00["COM"] = [s[0], s[1], COM[2]]

        float01 = copy.copy(float)
        float01["name"] = "float01"
        float01["translation"] = s = [0, shift_y/2., 0]
        float01["COM"] = [s[0], s[1], COM[2]]

        float02 = copy.copy(float)
        float02["name"] = "float02"
        float02["translation"] = s = [0, 0, 0]
        float02["COM"] = [s[0], s[1], COM[2]]

        float03 = copy.copy(float)
        float03["name"] = "float03"
        float03["translation"] = s = [0, -shift_y/2., 0]
        float03["COM"] = [s[0], s[1], COM[2]]

        inputfiles = [tank, water, float00, float01, float02, float03, absorber]

    elif "3floats" in suffix:

        '''
        y
        |--> x          

        00

        01

        02
        
        '''

        water["objfile"] = objfolder + "/water3floats.obj"

        shift_y = 300

        float00 = copy.copy(float)
        float00["name"] = "float00"
        float00["COM"] = [0, -shift_y, COM[2]]
        float00["translation"] = [0, -shift_y, 0]

        float01 = copy.copy(float)
        float01["name"] = "float01"
        float01["COM"] = [0, 0, COM[2]]

        float02 = copy.copy(float)
        float02["name"] = "float02"
        float02["COM"] = [0, shift_y, COM[2]]
        float02["translation"] = [0, shift_y, 0]

        inputfiles = [tank, water, float00, float01, float02, absorber]


    elif "2floats" in suffix:

        '''
        y
        |--> x          

        00

        01
        
        '''

        water["objfile"] = objfolder + "/water2floats.obj"

        shift_y = 300

        float00 = copy.copy(float)
        float00["name"] = "float00"
        float00["translation"] = s = [0, shift_y, 0]
        float00["COM"] = [s[0], s[1], COM[2]]

        float01 = copy.copy(float)
        float01["name"] = "float01"
        float01["COM"] = [0, 0, COM[2]]

        inputfiles = [tank, water, float00, float01, absorber]


    else:
        water["objfile"] = objfolder + "/water1float.obj"
        inputfiles = [tank, water, float, absorber]

    setting = {"max_dt": dt,
               "end_time_step": 10000000,
               "end_time": 100,
               "element": element,
               "ALE": ALE,
               "ALEPERIOD": ALEPERIOD}

    id = add_id(id)

    generate_input_files(inputfiles, setting, IO_dir, id)
elif "Tonegawa2024Akita" in SimulationCase:

    '''DOC_EXTRACT 2_1_0_input_generator

    ## Tonegawa2024 Akita

    ```sh
    python3.11 input_generator.py Tonegawa2024Akita -dt 0.2 -T 5 -H 2 -o /Volumes/home/BEM/Tonegawa2024Akita/ -s 3MW_MP30
    ```

    | 項目               | 重量 (ton)   | 主慣性モーメント(Ixx,Iyy,Izz) (kg-m²) |
    |-------------------|--------------|---------------------|
    |3MW_MP30x30        | 8437.5    | (0.67, 0.67, 0.995) x 1e10 |
    |3MW_MP15x15        | 13500     | (0.974, 0.974, 1.65) x 1e10 |
    |3MW_MP9x9          | 14580     | (1.04, 1.04, 1.75) x 1e10 |

    |10MW_MP30x30       | 13500     | (1.11,1.11,1.59) x 1e10 | 
    |10MW_MP15x15       | 21600     | (1.67,1.67,2.64) x 1e10 |
    |10MW_MP9x9         | 23328     | (1.77,1.77,2.83) x 1e10 |

    3MWの浮体のサイズは全て，WxLxH = 45x45x10 mで，喫水は7.5 m．
    10MWの浮体のサイズは全て，WxLxH = 60x60x17 mで，喫水は12 m．

    '''

    start = 0.

    a = H/2
    h = 80

    id = SimulationCase + "_H"+str(H).replace(".", "d")
    id += "_T"+str(T).replace(".", "d")

    if theta is not None:
        id += "_theta"+str(theta).replace(".", "d")

    water = {"name": "water", "type": "Fluid"}
    tank = {"name": "tank", "type": "RigidBody", "isFixed": True}
    float = {"name": "float", "type": "RigidBody", "velocity": "floating"}

    float_W = 45
    float_L = 45

    index = ""
    if "3MW_MP30" in suffix:
        draft = 7.5
        float_H = 10
        float["mass"] = m = 8437.5e3
        Ixx = Iyy = 0.67*1e10
        Izz = 0.995*1e10
        MP_W = MP_L = 30
        index = "3MW_MP30"
    elif "3MW_MP15" in suffix:
        draft = 7.5
        float_H = 10
        float["mass"] = m = 13500e3
        Ixx = Iyy = 0.974*1e10
        Izz = 1.65*1e10
        MP_W = MP_L = 15
        index = "3MW_MP15"
    elif "3MW_MP9" in suffix:
        draft = 7.5
        float_H = 10
        float["mass"] = m = 14580e3
        Ixx = Iyy = 1.04*1e10
        Izz = 1.75*1e10
        MP_W = MP_L = 9
        index = "3MW_MP9"
    elif "10MW_MP30" in suffix:
        draft = 12
        float_H = 17
        float["mass"] = m = 13500e3
        Ixx = Iyy = 1.11*1e10
        Izz = 1.59*1e10
        MP_W = MP_L = 30
        index = "10MW_MP30"
    elif "10MW_MP15" in suffix:
        draft = 12
        float_H = 17
        float["mass"] = m = 21600e3
        Ixx = Iyy = 1.67*1e10
        Izz = 2.64*1e10
        MP_W = MP_L = 15
        index = "10MW_MP15"
    elif "10MW_MP9" in suffix:
        draft = 12
        float_H = 17
        float["mass"] = m = 23328e3
        Ixx = Iyy = 1.77*1e10
        Izz = 2.83*1e10
        MP_W = MP_L = 9
        index = "10MW_MP9"
    else:
        print("Error: please specify the moon pool size")
        sys.exit()

    print("mass : ", m, "== mass based on the volume : ",
          1000*(float_W*float_L - MP_W*MP_L)*draft)

    MOI = [Ixx, Iyy, Izz]
    z_surface = 80
    z_floatinbody_bottom = z_surface - draft

    float["MOI"] = [Ixx, Iyy, Izz]
    float["COM"] = [0, 0, z_surface]

    objfolder = code_home_dir + "/cpp/obj/Tonegawa2024Akita"

    water["objfile"] = objfolder + "/water"+index+".obj"
    tank["objfile"] = objfolder + "/tank.obj"
    float["objfile"] = objfolder+"/float"+index+".obj"
    bottom_in_z = 0

    wave_theory = [a, T, h, bottom_in_z,
                   theta] if theta is not None else [a, T, h, bottom_in_z]

    H13 = H
    T13 = T
    random_wave_theory = [H13, T13, h, bottom_in_z]

    absorber = {"name": "absorber",
                "type": "Absorber",
                "isFixed": True,
                "objfile": objfolder+"/absorber.obj"}

    if "random_wave" in suffix:
        absorber["random_wave_theory"] = random_wave_theory
    else:
        absorber["wave_theory"] = wave_theory

    inputfiles = [tank, water, float, absorber]

    id = add_id(id)

    setting = {"max_dt": dt,
               "end_time_step": 10000000,
               "end_time": 100,
               "element": element,
               "ALE": ALE,
               "ALEPERIOD": ALEPERIOD}

    generate_input_files(inputfiles, setting, IO_dir, id)
elif "Ruehl2016" in SimulationCase:

    r'''DOC_EXTRACT 2_1_0_validation_FOSWEC

    ## Ruehl et al. (2016)
    
    査読付き学術論文ではないが，公開研究データセットと技術報告書

    Ruehl, K., Forbush, D., Lomonaco, P., Bosma, B., Simmons, A., Gunawan, B., Bacelli, G., & Michelen, C. (2016). Experimental Testing of a Floating Oscillating Surge Wave Energy Converter at Hinsdale Wave Research Laboratory. [Data set]. Marine and Hydrokinetic Data Repository. Sandia National Laboratories. https://doi.org/10.15473/1508364

    https://mhkdr.openei.org/submissions/310
    また，このデータセットは，Kelley Ruehl et al. (2020)のreferenceで引用されている．

    ```Mathematica
    n = 2000;
    timeData = 24*60*60 (import["/Users/tomoaki/Downloads/310/FOSWEC/data/WECSIM/inter/RegularWaveTuning1/Trial13/time.txt"][[;; n]] - import["/Users/tomoaki/Downloads/310/FOSWEC/data/WECSIM/inter/RegularWaveTuning1/Trial13/time.txt"][[1]]);
    dispData = import["/Users/tomoaki/Downloads/310/FOSWEC/data/WECSIM/inter/RegularWaveTuning1/Trial13/wmdisp15.txt"][[;; n]];
    data = Transpose[{timeData, dispData}];(*小数点10桁、指数なしで表示形式に変換*)
    cppList = StringJoin["{", StringRiffle[Map[StringJoin["{", StringRiffle[Map[ToString[NumberForm[#, {Infinity, 10}, ExponentFunction -> (Function[exp, Null])]] &, #], ","], "}"] &, data], ","], "}"];(*ファイルに書き出し*)
    Export["/Users/tomoaki/Downloads/data_for_cpp.txt", cppList, "Text"]
    ```
    
    Directional Wave Basisは
    Length: 48.8 m
    Width: 26.5 m
    Max depth: 1.37 m

    Wavemaker
    Type: Piston-type
    Waceboards: 29 boards, 2.0 m high
    Period range: 0.5 to 10 s
    
    項目	内容
    装置名	FOSWEC（Floating Oscillating Surge Wave Energy Converter）
    目的	モーション制御可能な2-DOF型WEC（Wave Energy Converter）の応答特性・発電性能評価
    試験場所	OSU Hinsdale Wave Research Laboratory
    模型スケール	1:33 スケールモデル
    自由動揺の自由度	Surge（前後）、Pitch（ピッチ）の2自由度
    波の条件	単一周波数（monochromatic）および擬似ランダム（pseudo-random）波を使用
    計測項目	本体動揺（surge, pitch）、アーム角、PTO（Power Take Off）からの出力電力、波高など
    PTOの形式	モータと可変抵抗による受動的な減衰器（トルク制御可能）
    主な結果	共振周波数において最大応答を示し、周波数応答関数（RAO）を明確に評価。動揺と発電出力の相関性を検証。RAOの実験値と数値モデルとの比較も可能。

    <img src="./FOSWEC_testReport_table38.png" alt="FOSWEC_testReport_table38.png" width="600px">
    <img src="./FOSWEC_testReport_table41.png" alt="FOSWEC_testReport_table41.png" width="600px">

    FOSWEC/
    ├── data/
    │   ├── WECSIM/    ← Phase 1 データ
    │   │   ├── logs/     ← 試験ログ（Appendix D）
    │   │   └── inter/    ← 中間処理済みデータ（個別実験別フォルダ）
    │   └── WECSIM2/   ← Phase 2 データ
    │       ├── logs/     ← 試験ログ（Appendix E）
    │       └── inter/    ← 中間処理済みデータ（個別実験別フォルダ）
    ├── doc/            ← 説明資料
    └── geom/           ← 幾何形状データ

    ### 水槽における計測位置

    wg3.txtなどにコメントとして書いてある

    WG3: (x,y) = (9.480,0)
    wg6: (x,y) = (16.778, -0.031)

    水深: h = 1.36 m
                            
    '''

    # IDの生成
    id = SimulationCase
    id = add_id(id)
    max_dt = dt

    objfolder = code_home_dir + "/cpp/obj/Ruehl2016/"

    water = {"name": "water", "type": "Fluid", "objfile": objfolder + "/water.obj"}

    if meshname != "":
        water["objfile"] = objfolder + "/" + meshname

    vertices, triangles = read_obj(water["objfile"])
    water["vertices"] = vertices
    water["triangles"] = triangles

    tank = {"name": "tank", "type": "RigidBody", "isFixed": True, "objfile": objfolder + "/tank.obj"}

    h = 1.36
    # wavemaker = {"name": "wavemaker", "type": "RigidBody", "velocity": ["Goring1979", 3., 0.1*h, h], "objfile": f"{objfolder}/wavemaker.obj"}
    trial = suffix
    wavemaker = {"name": "wavemaker", "type": "RigidBody", "velocity": ["file", "./study_benchmark/Ruehl2016/"+ trial +".csv"], "objfile": f"{objfolder}/wavemaker.obj"}

    absorber = {"name": "absorber", "type": "Absorber", "isFixed": True, "objfile": objfolder + "/absorber.obj"}

    if T is not None:
        absorber["wave_theory"] = [0, T, h, 0]
    elif L is not None:
        absorber["wave_theory_L"] = [0, L, h, 0]
    else:
        absorber["wave_theory_L"] = [0, 1., h, 0]

    gauges = []
    gauges.append({"name": "wg3", "type": "wave gauge", "position": [9.48, 0., 1.36+0.5, 9.48, 0., 1.36-0.5]})
    gauges.append({"name": "wg6", "type": "wave gauge", "position": [16.778, -0.031, 2., 16.778, -0.031, 0.5]})

    inputfiles = [tank, wavemaker, water] + gauges
    
    setting = {"max_dt": max_dt,
               "end_time_step": 100000,
               "end_time": 35,
               "element": element,
               "ALE": ALE,
               "ALEPERIOD": ALEPERIOD}

    generate_input_files(inputfiles, setting, IO_dir, id)

else :
    # シミュレーションケースの指定がない場合はエラー
    # SimulationCase is
    print("SimulationCase", SimulationCase)
    print("Error: please specify the simulation case")
    sys.exit()
# シミュレーションケースのモジュールを動的にインポート
# case_module = importlib.import_module(f"cases.{args.case}")
# case_module.generate_input_files(args)