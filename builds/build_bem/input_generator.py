'''DOC_EXTRACT 2_0_0_input_generator

# Input Generator

This file is used to generate the input files for the BEM-MEL.

'''
import copy, platform, json, math, os, sys
from math import pi
from os.path import expanduser

sys.path.append("..")
from IOgenerator import generate_input_files

input_dir = "./input_files/"
sys_home_dir = expanduser("~")
current_dir = os.path.dirname(os.path.abspath(__file__))
code_home_dir = os.path.join(current_dir, '../../../../code')

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

match SimulationCase:

    case "fish_without_free_surface":

        start = 0.

        id = SimulationCase

        objfolder = code_home_dir + "/cpp/obj/fish"

        waterA = {"name": "waterA",
                 "type": "Fluid",
                 "reverseNormal" : True,
                 "objfile": objfolder + "/bodyA50.obj"}

        waterB = {"name": "waterB",
                 "type": "Fluid",
                 "reverseNormal" : True,
                 "objfile": objfolder + "/bodyB50.obj"}

        waterC = {"name": "waterC",
                 "type": "Fluid",
                 "reverseNormal" : True,
                 "objfile": objfolder + "/bodyC50.obj"}

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

        inputfiles = [waterA, waterB, waterC, bodyA, bodyB, bodyC]

        setting = {"max_dt": 0.01,
                   "end_time_step": 10000,
                   "end_time": 9}

        generate_input_files(inputfiles, setting, IO_dir, id)
    case "fish":

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
    case "testALE":

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
    case "Ren2015":

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
        H = 0.06
        a = H/2  # Ren 2015 used H=[0.1(a=0.05), 0.03(a=0.06), 0.04(a=0.02)]
        h = 0.4

        id0 = ""
        # id0 = "_no"
        # id0 = "_multiple"

        wavemaker_type = "piston"
        # wavemaker_type = "potential"

        id = SimulationCase + id0 + "_H"+str(H).replace(".", "d")
        id += "_T"+str(T).replace(".", "d")
        id += "_" + wavemaker_type + "_2"

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
    case "Hadzic2005":
        
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

        multiple_case = False

        objfolder = code_home_dir + "/cpp/obj/Hadzic2005_24floats" if multiple_case else code_home_dir + "/cpp/obj/Hadzic2005"

        water = {"name": "water", 
                 "type": "Fluid",
                 "objfile": objfolder + "/water1000meshlab.obj" if multiple_case else objfolder + "/water1000meshlab3.obj"}
        
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
                 "velocity": "floating"
                 }
        
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
        float["MOI"] = [Ixx*10**10,0.01*0.01*20,Izz*10**10]
        inputfiles = [tank, wavemaker, water, float]

        setting = {"max_dt": 0.05,
                   "end_time_step": 10000,
                   "end_time": 9}

        if multiple_case:
            id = SimulationCase + "_multiple"
        else:
            id = SimulationCase

        generate_input_files(inputfiles, setting, IO_dir, id)
    case "Kramer2021":

        '''DOC_EXTRACT 2_1_1_validation_Kramer2021
        
        This case is for the validation of the floating body motion analysis using the BEM-MEL.

        <img src="schematic_Kramer2021.png" width="400px" />

        The floating body is a sphere with the diameter of 0.3 m.
        The mass of the floating body is 7.056 kg.
        The moment of inertia of the floating body is set to be almost infinite to ignore the effect of the rotation.

        The sphere is dropped from the height of 0.03 m above the water surface.
        '''
        start_t = max_dt = 0.02
        D = 300/1000
        H0 = D*0.1
        id = "H0"+str(H0).replace(".", "d") + "_small"
        # id = "H0"+str(H0).replace(".", "d")

        input_dir += SimulationCase + "_" + id
        os.makedirs(input_dir, exist_ok=True)
        output_dir = sys_home_dir + "/BEM/" + SimulationCase + "_" + id
        os.makedirs(output_dir, exist_ok=True)

        water = {"name": "water", "type": "Fluid"}
        tank = {"name": "tank", "type": "RigidBody", "isFixed": True}

        float = {"name": "float",
                 "type": "RigidBody",
                 "velocity": ["floating", start_t]}

        float["mass"] = m = 7.056
        # float["reverseNormal"] = True
        z_surface = 900/1000
        float["COM"] = [0., 0., z_surface + 0.1*D]  # 今回は重要ではない
        float["radius_of_gyration"] = [10**10, 10**10, 10**10]
        float["MOI"] = [m*math.pow(float["radius_of_gyration"][0], 2),
                        m*math.pow(float["radius_of_gyration"][1], 2),
                        m*math.pow(float["radius_of_gyration"][2], 2)]
        # float["translate"] = [0., 0., 0.1*D + 900/1000]

        objfolder = code_home_dir + "/cpp/obj/" + SimulationCase + "_" + id

        if id == "H00d03":
            water["objfile"] = objfolder + "/water200_mod.obj"
        else:
            water["objfile"] = objfolder + "/water300_mod.obj"

        tank["objfile"] = objfolder + "/tank10.obj"
        float["objfile"] = objfolder + "/sphere.obj"

        inputfiles = [tank, water, float]

        rho = 998.2
        g = 9.82
        setting = {"WATER_DENSITY": rho,
                   "GRAVITY": g,
                   "max_dt": max_dt,
                   "end_time_step": 10000,
                   "end_time": 4}

        id = SimulationCase + "_" + id
        generate_input_files(inputfiles, setting, IO_dir, id)
    case "simple_barge":

        input_dir += SimulationCase
        os.makedirs(input_dir, exist_ok=True)
        output_dir = sys_home_dir + "/BEM/" + SimulationCase
        os.makedirs(output_dir, exist_ok=True)

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
                         "type": "SoftBody",
                         "isFixed": True,
                         "velocity": ["linear_traveling_wave", start, a, T, h, z_surface]}

        if False:
            floatingbody = {"name": "floatingbody",
                            "type": "RigidBody",
                            "velocity": ["sin", 0, a, T]}
        else:
            floatingbody = {"name": "floatingbody",
                            "type": "RigidBody",
                            "velocity": "floating"}

        A = 2450.00
        floatingbody["mass"] = m = (rho*g*7.5*A)/g
        floatingbody["COM"] = [100., 50., 75.]
        floatingbody["radius_of_gyration"] = [20., 20., 20.]
        floatingbody["MOI"] = [m*math.pow(floatingbody["radius_of_gyration"][0], 2),
                               m*math.pow(floatingbody["radius_of_gyration"][1], 2),
                               m*math.pow(floatingbody["radius_of_gyration"][2], 2)]

        objfolder = code_home_dir + "/cpp/obj/tsukada2022_re"
        water["objfile"] = objfolder + "/water100_mod.obj"
        wavemaker["objfile"] = objfolder + "/wavemaker200.obj"
        tank["objfile"] = objfolder + "/wavetank10.obj"
        floatingbody["objfile"] = objfolder+"/floatingbody100.obj"

        inputfiles = [tank, wavemaker, water, floatingbody]

        setting = {"max_dt": 0.25,
                   "end_time_step": 10000,
                   "end_time": 10000}

        id = SimulationCase
        generate_input_files(inputfiles, setting, IO_dir, id)
    case "moon_pool":

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

            floatingbody = {"name": "floatingbody",
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
    case "two_floatingbodies":

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

        floatingbody_a = {"name": "floatingbody_a",
                          "COM": [300., 150., 150-78./2]}

        floatingbody_b = {"name": "floatingbody_b",
                          "COM": [600., 150., 150-78./2]}

        floatingbody_c = {"name": "floatingbody_c",
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
    case "Gu2018Float1_offset_neg0d074":

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
    case "Tonegawa":

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

        floatingbody = {"name": "floatingbody",
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
    case "three_floatingbodies":

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

        floatingbody_a = {"name": "floatingbody_a",
                          "COM": [300., 150., 150-78./2]}

        floatingbody_b = {"name": "floatingbody_b",
                          "COM": [600., 150., 150-78./2]}

        floatingbody_c = {"name": "floatingbody_c",
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
    case "Retzler2000simple":

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
                     "isFixed": True,
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
    case "Retzler2000":

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
                     "isFixed": True,
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
    case "three_200":

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
