import os
from IOgenerator import generate_input_files, read_obj

def generate_input_files(args):
    objfolder = os.path.join(os.path.dirname(__file__), "../../../cpp/obj/Cheng2018")

    water = {"name": "water", "type": "Fluid"}
    tank = {"name": "tank", "type": "RigidBody", "isFixed": True}
    absorber = {"name": "absorber", "type": "Absorber", "isFixed": True}

    start = 0.
    T = 1.2
    H = 0.04
    a = H / 2
    h = 0.4

    if "meshA" in args.case:
        id0 = "meshA"
    elif "meshB" in args.case:
        id0 = "meshB"
    elif "meshC" in args.case:
        id0 = "meshC"
    else:
        id0 = "meshA"

    wavemaker_type = "piston"
    id = args.case
    id += "_H" + str(H).replace(".", "d")
    id += "_T" + str(T).replace(".", "d")
    id += "_" + wavemaker_type
    id += "_with_float"

    gauges = []
    dx = 0.5
    for i in range(10):
        gauges.append({"name": "gauge" + str(i), "type": "wave gauge", "position": [0.2 + dx * i, 0, 0.6, 0.2 + dx * i, 0, 0.1]})

    z_surface = 0.4
    phase_shift = -pi / 2.
    if "free_decay" in id:
        wavemaker = {"name": "wavemaker", "type": "RigidBody", "velocity": ["piston", start, 0, T, h, 1, 0, 0]}
    elif "piston" in wavemaker_type:
        wavemaker = {"name": "wavemaker", "type": "RigidBody", "velocity": ["piston", start, a, T, h, 1, 0, 0]}
    elif "flap" in wavemaker_type:
        wavemaker = {"name": "wavemaker", "type": "RigidBody", "velocity": ["flap", start, a, T, h, h, 0, 1, 0]}
    else:
        wavemaker = {"name": "wavemaker", "type": "SoftBody", "isFixed": True, "velocity": ["velocity", start, a, T, h, z_surface, phase_shift]}

    Lx = 0.3
    Lz = 0.2
    d = 0.1
    Ly = 0.42

    float = {"name": "float", "type": "RigidBody", "velocity": "floating", "isFixed": [False, True, False, True, False, True]}
    float["mass"] = m = 500 * Lx * Lz * Ly
    Ixx = 1. / 12. * m * (Ly * Ly + Lz * Lz)
    Iyy = 0.1967978007
    Izz = 1. / 12. * m * (Lx * Lx + Ly * Ly)
    z_surface = 0.4
    float["COM"] = [3.35, 0., z_surface - d]
    float["MOI"] = [10 ** 10 * Ixx, Iyy, 10 ** 10 * Izz]

    if "free_decay" in id:
        water["objfile"] = os.path.join(objfolder, "water_inclined_mod.obj")
        inputfiles = [tank, wavemaker, water, float, absorber]
        float["objfile"] = os.path.join(objfolder, "float_inclined.obj")
        wavemaker["velocity"] = ["velocity", start, 0., T, h, z_surface, phase_shift]
    elif "without_float" in id:
        water["objfile"] = os.path.join(objfolder, "water_without_float9_mod.obj")
        inputfiles = [tank, wavemaker, water, absorber]
        float["objfile"] = os.path.join(objfolder, "float7.obj")
    else:
        if id0 != "":
            water["objfile"] = os.path.join(objfolder, "water8_" + id0 + ".obj")
        else:
            water["objfile"] = os.path.join(objfolder, "water8_modmodmod.obj")
        inputfiles = [tank, wavemaker, water, float, absorber]
        float["objfile"] = os.path.join(objfolder, "float7.obj")

    wavemaker["objfile"] = os.path.join(objfolder, "wavemaker8.obj")
    tank["objfile"] = os.path.join(objfolder, "tank4.obj")
    absorber["objfile"] = os.path.join(objfolder, "absorber.obj")
    inputfiles += gauges

    setting = {"max_dt": 0.02, "end_time_step": 10000, "end_time": 20}

    generate_input_files(inputfiles, setting, lambda id: (f"./input_files/{id}", f"{args.outputdir}/{id}"), id)
