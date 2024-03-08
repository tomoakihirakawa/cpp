'''DOC_EXTRACT 2_0_0_input_generator

# Input Generator

This file is used to generate the input files for the BEM-MEL.

'''

import copy
import platform
import json
import math
import os
import sys
from math import pi
from os.path import expanduser

sys.path.append("..")

input_dir = "./input_files/"
sys_home_dir = expanduser("~")
current_dir = os.path.dirname(os.path.abspath(__file__))
code_home_dir = os.path.join(current_dir, '../../../../code')
from IOgenerator import generate_input_files, read_obj, Norm3d, Differece3d, Add3d


def IO_dir(id):
    input_dir = "./input_files/" + id
    os.makedirs(input_dir, exist_ok=True)
    output_dir = sys_home_dir + "/BEM/" + id
    # output_dir = "/Volumes/home/BEM/benchmark20231206/" + id
    os.makedirs(output_dir, exist_ok=True)
    return input_dir, output_dir


rho = 1000.
g = 9.81

# ---------------------------------------------------------------------------- #

SimulationCase = "Hadzic2005"
if len(sys.argv) > 1:
    # If a value is passed, use it
    SimulationCase = sys.argv[1]
    print(f"Value passed from command line: {SimulationCase}")
else:
    # Default action if no value is passed
    print("No value passed. Executing default action.")

# ---------------------------------------------------------------------------- #
if "Tanizawa1996" in SimulationCase:

    '''DOC_EXTRACT 2_1_0_validation_Tanizawa1996

    <img src="schematic_float_Tanizawa1996.png" width="400px" />

    \cite{Tanizawa1996}には，BEM-MELを使う際の消波方法が紹介されている．
    この理論は簡単で，$`\frac{D\phi}{Dt}`$と$`\frac{D\bf x}{Dt}`$に減衰項を追加するだけである．
    興味深いのは，この項は，単なる消波だけでなく，ある領域の表面変位と速度ポテンシャルを理想的な値へと保つことができるため，造波装置側に設置するのも有効であることである．

    浮体の左右には整流板が設置されている．

    水深は4.5mで一定．

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

    L = 1.8
    H = 0.05
    a = H/2
    h = 5.
    z_surface = h

    # wavemaker_type = "piston"
    wavemaker_type = "potential"
    # wavemaker_type = "flap"

    id = SimulationCase

    meshname = "water0d09"
    # meshname = "water_mod_0d05"
    # meshname = "multiple_water_mod_coarse"
    id += "_H" + str(H).replace(".", "d")
    id += "_L" + str(L).replace(".", "d")
    id += "_" + wavemaker_type
    id+= "_" +  meshname

    water = {"name": "water", "type": "Fluid"}
    tank = {"name": "tank", "type": "RigidBody", "isFixed": True}
    absorber = {"name": "absorber", "type": "Absorber", "isFixed": True}

    gauges = []
    dx = 0.5
    for i in range(24):
        gauges.append({"name": "gauge"+str(i),
                       "type": "wave gauge",
                       "position": [0.2 + dx*i, 0, h + 2, 0.2 + dx*i, 0, h - 2]})

    phase_shift = -pi/2.
    # if wavemaker_type contains "piston":

    if "piston" in wavemaker_type:
        wavemaker = {"name": "wavemaker",
                    "type": "RigidBody",
                    "velocity": ["piston_using_wave_length", start, a, L, h, 1, 0, 0]}
    elif "potential" in wavemaker_type:
        wavemaker = {"name": "wavemaker",
                    "type": "SoftBody",
                    "isFixed": True,
                    "velocity": ["velocity__using_wave_length", start, a, L, h, h]}

    Lx = 0.74
    Lz = 0.415
    d = 0.25
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
                        "isFixed": [False,True,False,True,False,True]}
                m = 184.3
                float["mass"] = m
                Ixx = 1./12.*m*(Ly*Ly+Lz*Lz)
                Iyy = m*math.pow(0.266, 2)
                Izz = 1./12.*m*(Lx*Lx+Ly*Ly)
                float["COM"] = [12 - 1.5*(i-1),-2 + 2 *(j-1) ,z_surface - d + 0.22]
                float["spring"] = [float["COM"][0], float["COM"][1], float["COM"][2], 51., 1000., 0.]
                float["MOI"] = [1000*Iyy, Iyy, 1000*Iyy]
                float["objfile"] = objfolder+"/float"+str(i+8*(j-1))+".obj"
                floats.append(float)
        water["objfile"] = objfolder + "/water_mod_coarse.obj"
    else:
        objfolder = code_home_dir + "/cpp/obj/Tanizawa1996_2d"
        m = 184.3
        Ixx = 1./12.*m*(Ly*Ly+Lz*Lz)
        K = 0.266
        Iyy = m*math.pow(K, 2) # 0.266 is the radius of inertia
        Izz = 1./12.*m*(Lx*Lx+Ly*Ly)
        float = {"name": "float",
                "type": "RigidBody",
                "velocity": "floating",
                "damping": [1000, 1000, 1000, 1000, 1000, 1000, 0., 11.],
                "isFixed": [False,True,False,True,False,True],
                "mass" : m}

        # 密度1.18g/cm^3と質量から計算，{Ixx,Iyy,Izz}={0.3034752725,0.1967978007,0.3692731657}
        # 板厚0.08mmと質量から計算した結果，{Ixx,Iyy,Izz}={0.3320104413,0.2204711483,0.4018598521}

        # float["COM"] = [7.5, 0., z_surface - d + 0.22]
        float["COM"] = [12., 0., z_surface - d + 0.22]
        float["spring"] = [float["COM"][0], float["COM"][1], float["COM"][2], 51., 1000., 0.]#実際は[51., 0., 0.]だがy方向にずれが生じるためこのようにした
        print("COM ", float["COM"], " mass ", float["COM"], " Ly ", Ly)
        float["MOI"] = [1000*Iyy, Iyy, 1000*Iyy]
        float["objfile"] = objfolder+"/float.obj"
        floats.append(float)

        # water["objfile"] = objfolder + "/water_mod1.obj"
        water["objfile"] = objfolder + "/" + meshname + ".obj"

    # water["objfile"] = objfolder + "/water_mod_coarse.obj"
    
    wavemaker["objfile"] = objfolder + "/wavemaker.obj"

    tank["objfile"] = objfolder + "/tank.obj"
    absorber["objfile"] = objfolder+"/absorber.obj"

    inputfiles = [tank, wavemaker, water, absorber]
    inputfiles += gauges
    inputfiles += floats

    setting = {"max_dt": 1/40.,
               "end_time_step": 20000,
               "end_time": 30}

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
             "isFixed": [False,True,False,True,False,True]}

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
        wavemaker["velocity"] = ["velocity", start, 0., T, h, z_surface, phase_shift]
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

    '''DOC_EXTRACT 2_1_0_validation_Ren2015

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
    id += "_" + wavemaker_type

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
    H = 0.2
    d = 0.1
    # 500*L*W*H (浮体全質量) = 1000*L*W*d (排除される水の質量)
    # d = H/2 = 0.1
    W = 0.42

    float = {"name": "float",
             "type": "RigidBody",
             # "isFixed": True,
             # "output": "json"}
             # "mooring": ["simple_mooring", 4.6, W/2, 0., 4.6, W/2., 0.3, 1.],
             #  "mooring": ["mooring1", # 1
             #              "simple_mooring", # 2
             #              4.6, W/2, 0.,   # x,y,z body side
             #              4.6, W/2., 0.3, # x,y,z
             #              10,     #div 9
             #              10**7., # stiffness 10
             #              0.9,    # damping 11
             #              348.5   # weight 12
             #              ],
             "velocity": "floating"}

    float["mass"] = m = 500*L*H*W
    Ixx = 1./12.*m*(W*W+H*H)
    Iyy = 1./12.*m*(L*L+H*H)
    Iyy = 0.1967978007  # MathematicaでrhoPMMA=1180 kg/m^3 として計算した結果
    Izz = 1./12.*m*(L*L+W*W)
    MOI = [Ixx, Iyy, Izz]
    z_surface = 0.4
    z_floatinbody_bottom = z_surface - d
    # float["COM"] = [2+L/2., W/2, z_surface]
    float["COM"] = [4.6+L/2, W/2, z_surface]
    print("COM ", float["COM"], " mass ", float["COM"], " W ", W)
    float["MOI"] = [10**10*Ixx, Iyy, 10**10*Izz]

    # if id contains "multiple":
    if "multiple" in id:
        float["COM"] = [4.3+L/2, W/2, z_surface]
        objfolder = code_home_dir + "/cpp/obj/Ren2015_multiple"

        # Initialize the object files
        water["objfile"] = f"{objfolder}/water400.obj"
        wavemaker["objfile"] = f"{objfolder}/wavemaker20.obj"
        tank["objfile"] = f"{objfolder}/tank50.obj"
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
        tank["objfile"] = objfolder + "/tank10.obj"
        inputfiles = [tank, wavemaker, water]
    else:
        float["COM"] = [4.6, W/2, z_surface]
        objfolder = code_home_dir + "/cpp/obj/Ren2015"
        # water["objfile"] = objfolder + "/water400meshlab.obj"
        water["objfile"] = objfolder + "/water18_mod.obj"
        wavemaker["objfile"] = objfolder + "/wavemaker20.obj"
        tank["objfile"] = objfolder + "/tank10.obj"
        float["objfile"] = objfolder+"/float10.obj"
        inputfiles = [tank, wavemaker, water, float]

    setting = {"max_dt": 0.01,
               "end_time_step": 10000,
               "end_time": 9}

    generate_input_files(inputfiles, setting, IO_dir, id)
elif "Hadzic2005" in SimulationCase:

    '''DOC_EXTRACT 2_1_0_validation_Hadzic2005                

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

    '''DOC_EXTRACT 2_1_1_validation_Kramer2021

    This case is for the validation of the floating body motion analysis using the BEM-MEL.

    <img src="schematic_Kramer2021.png" width="400px" />

    The floating body is a sphere with the diameter of 0.3 m.
    The mass of the floating body is 7.056 kg.
    The moment of inertia of the floating body is set to be almost infinite to ignore the effect of the rotation.

    The sphere is dropped from the height of 0.03 m above the water surface.
    '''

    max_dt = 0.001
    D = 300/1000
    H0 = D*0.3
    small = True
    id = "H0"+str(H0).replace(".", "d")
    if small:
        id += "_small"
    # id = "H0"+str(H0).replace(".", "d")

    input_dir += SimulationCase + "_" + id
    os.makedirs(input_dir, exist_ok=True)
    output_dir = sys_home_dir + "/BEM/" + SimulationCase + "_" + id
    os.makedirs(output_dir, exist_ok=True)

    water = {"name": "water", "type": "Fluid"}
    tank = {"name": "tank", "type": "RigidBody", "isFixed": True}

    start_t = 0.02
    float = {"name": "float",
             "type": "RigidBody",
             "velocity": ["floating", start_t]}

    float["mass"] = m = 7.056
    # float["reverseNormal"] = True
    z_surface = 900/1000
    float["COM"] = [0., 0., z_surface + H0]  # 今回は重要ではない
    float["radius_of_gyration"] = [10**10, 10**10, 10**10]
    float["MOI"] = [m*math.pow(float["radius_of_gyration"][0], 2),
                    m*math.pow(float["radius_of_gyration"][1], 2),
                    m*math.pow(float["radius_of_gyration"][2], 2)]
    # float["translate"] = [0., 0., 0.1*D + 900/1000]

    mesh = "meshB"
    id += "_" + mesh
    if "H00d03" in id:
        objfolder = code_home_dir + "/cpp/obj/" + SimulationCase + "_01D"
        if small:
            objfolder += "_small"
        water["objfile"] = objfolder + "/water5_" + mesh + ".obj"
    elif "H00d09" in id:
        objfolder = code_home_dir + "/cpp/obj/" + SimulationCase + "_03D"
        if small:
            objfolder += "_small"
        water["objfile"] = objfolder + "/water4_" + mesh + ".obj"

    tank["objfile"] = objfolder + "/tank4.obj"
    float["objfile"] = objfolder + "/sphere.obj"

    gauges = []
    dx = 0.6
    for i in range(10):
        gauges.append({"name": "gauge"+str(i), "type": "wave gauge",
                      "position": [0. + dx*i, 0, 1.1, 0. + dx*i, 0, 1.1-0.6]})

    inputfiles = [tank, water, float]
    inputfiles += gauges

    rho = 998.2
    g = 9.82

    setting = {"WATER_DENSITY": rho,
               "GRAVITY": g,
               "max_dt": max_dt,
               "end_time_step": 10000,
               "end_time": 4,
               "ALE_period_in_step": 1}

    id = SimulationCase + "_" + id
    id += "_ALE" + str(setting["ALE_period_in_step"]) + "_0d5"
    generate_input_files(inputfiles, setting, IO_dir, id)
elif "Palm2016" in SimulationCase:

    '''DOC_EXTRACT 2_1_0_validation_Palm2016

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
    a = 0.02
    H = 2*a
    T = 1.
    h = 0.9
    z_surface = 0.9

    SimulationCase += "_H" + str(H).replace(".", "d")
    SimulationCase += "_T" + str(T).replace(".", "d")

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
    X_fair_leadA = [r * math.cos(q) + COM[0], r*math.sin(q) + COM[1], z_surface]
    X_anchorA = Add3d(X_fair_leadA, [horizontal_length*math.cos(q), horizontal_length*math.sin(q), -z_surface])

    q = (i+1)*2*pi/3
    X_fair_leadB = [r * math.cos(q) + COM[0], r*math.sin(q) + COM[1], z_surface]
    X_anchorB = Add3d(X_fair_leadB, [horizontal_length*math.cos(q), horizontal_length*math.sin(q), -z_surface])

    q = (i+2)*2*pi/3
    X_fair_leadC = [r * math.cos(q) + COM[0], r*math.sin(q) + COM[1], z_surface]
    X_anchorC = Add3d(X_fair_leadC, [horizontal_length*math.cos(q), horizontal_length*math.sin(q), -z_surface])

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
             "MOI": [Ixx, Ixx, 10.**10]}

    id = SimulationCase

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
    water["objfile"] = objfolder + "/water_mod.obj"
    wavemaker["objfile"] = objfolder + "/wavemaker_mod.obj"
    tank["objfile"] = objfolder + "/tank_mod.obj"
    float["objfile"] = objfolder+"/float_mod.obj"

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
elif "Liang2022" in SimulationCase:

    '''DOC_EXTRACT 2_1_0_validation_Liang2022

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
        r = math.sqrt(1.567**2 - float_bottom**2) # horizontal_length
        X_fair_leadA = [float["COM"][0] - Lx/2, float["COM"][1] - Ly/2, float_bottom]
        X_anchorA = [float["COM"][0] - Lx/2 + r, float["COM"][1] - Ly/2, 0]
        X_fair_leadB = [float["COM"][0] - Lx/2, float["COM"][1] + Ly/2, float_bottom]
        X_anchorB = [float["COM"][0] - Lx/2 + r, float["COM"][1] + Ly/2, 0]
        X_fair_leadC = [float["COM"][0] + Lx/2, float["COM"][1] - Ly/2, float_bottom]
        X_anchorC = [float["COM"][0] + Lx/2 - r, float["COM"][1] - Ly/2, 0]
        X_fair_leadD = [float["COM"][0] + Lx/2, float["COM"][1] + Ly/2, float_bottom]
        X_anchorD = [float["COM"][0] + Lx/2 - r, float["COM"][1] + Ly/2, 0]
    elif "typeB" in SimulationCase:
        r = math.sqrt(1.125**2 - float_bottom**2)
        X_fair_leadA = [float["COM"][0] - Lx/2, float["COM"][1] - Ly/2, float_bottom]
        X_anchorA = [float["COM"][0]  - Lx/2 - r, float["COM"][1] - Ly/2, 0]
        X_fair_leadB = [float["COM"][0] - Lx/2, float["COM"][1] + Ly/2, float_bottom]
        X_anchorB = [float["COM"][0]  - Lx/2 - r, float["COM"][1] + Ly/2, 0]
        X_fair_leadC = [float["COM"][0] + Lx/2, float["COM"][1] - Ly/2, float_bottom]
        X_anchorC = [float["COM"][0]  + Lx/2 + r, float["COM"][1] - Ly/2, 0]
        X_fair_leadD = [float["COM"][0] + Lx/2, float["COM"][1] + Ly/2, float_bottom]
        X_anchorD = [float["COM"][0]  + Lx/2 + r, float["COM"][1] + Ly/2, 0]
    elif "typeC" in SimulationCase:
        r = math.sqrt(0.809**2 - float_bottom**2)
        X_fair_leadA = [float["COM"][0] - Lx/2, float["COM"][1] - Ly/2, float_bottom]
        X_anchorA = [float["COM"][0] - Lx/2 - r, float["COM"][1] - Ly/2, 0]
        X_fair_leadB = [float["COM"][0] - Lx/2, float["COM"][1] + Ly/2, float_bottom]
        X_anchorB = [float["COM"][0] - Lx/2- r, float["COM"][1] + Ly/2, 0]
        X_fair_leadC = [float["COM"][0] + Lx/2, float["COM"][1] - Ly/2, float_bottom]
        X_anchorC = [float["COM"][0] + Lx/2 + r, float["COM"][1] - Ly/2, 0]
        X_fair_leadD = [float["COM"][0] + Lx/2, float["COM"][1] + Ly/2, float_bottom]
        X_anchorD = [float["COM"][0] + Lx/2+ r, float["COM"][1] + Ly/2, 0]

    stiffness = 2.36*10**3  # ! [N/m]
    damp = .4  # ! [N/(m/s^2)]
    density = 0.177  # ! [kg/m]

    total_lengthA = Norm3d(Differece3d(X_anchorA,X_fair_leadA))
    total_lengthB = Norm3d(Differece3d(X_anchorB,X_fair_leadB))
    total_lengthC = Norm3d(Differece3d(X_anchorC,X_fair_leadC))
    total_lengthD = Norm3d(Differece3d(X_anchorD,X_fair_leadD))

    n_points = 30
    diam = 1*0.1 # 1 cm

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
    id += '_with_mooring'

    water = {"name": "water", "type": "Fluid"}

    tank = {"name": "tank", "type": "RigidBody", "isFixed": True}

    start = 0.
    a = 2.
    T = 7.
    h = 80
    z_surface = 80

    if False:
        wavemaker = {"name": "wavemaker",
                     "type": "RigidBody",
                     "isFixed": True}
    else:
        wavemaker = {"name": "wavemaker",
                     "type": "RigidBody",
                     "velocity": ["piston", start, a, T, h, 1, 0, 0]}

    # 係留索の設定

    X_anchorA = [50., 50., 0.]
    X_fair_leadA = [75., 50., 71.]

    X_anchorB = [150., 50., 0.]
    X_fair_leadB = [125., 50, 71.]

    stiffness = 14 * 10**8  # ! [N/m]
    damp = .5  # ! [N/(m/s^2)]
    density = 100.  # ! [kg/m]
    dt = 0.01  # ! [s]
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

    setting = {"max_dt": 0.1,
               "end_time_step": 100000,
               "end_time": 100}

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
                     "velocity": ["linear_traveling_wave", start, a, T, h, z_surface]}  # "isFixed": True}

        floatingbody = {"name": "float",
                        "type": "RigidBody",
                        "velocity": "floating"}  # "velocity": ["sin", 0, a, T]}

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

        setting = {"max_dt": 0.2,
                   "end_time_step": 1000,
                   "end_time": 100}

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
elif "Tonegawa" in SimulationCase:
    input_dir += SimulationCase
    os.makedirs(input_dir, exist_ok=True)
    output_dir = sys_home_dir + "/BEM/2022Tonegawa3"
    os.makedirs(output_dir, exist_ok=True)

    water = {"name": "water",
             "type": "Fluid"
             }

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
    setting = {"max_dt": 0.002,
               "end_time_step": 1000,
               "end_time": 0.4}

    id = SimulationCase
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
