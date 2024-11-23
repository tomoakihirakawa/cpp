import os
from IOgenerator import generate_input_files, read_obj

def generate_input_files(args):
    objfolder = os.path.join(os.path.dirname(__file__), "../../../cpp/obj/Tanizawa1996_2d")

    water = {"name": "water", "type": "Fluid"}
    tank = {"name": "tank", "type": "RigidBody", "isFixed": True}
    absorber = {"name": "absorber", "type": "Absorber", "isFixed": True}

    start = 0.
    L = 1.8
    if args.H is None:
        H = 0.05 * 2
    else:
        H = args.H

    a = H / 2
    h = 5.
    z_surface = h

    id = args.case

    if args.wavemaker == "":
        wavemaker_type = "potential"
    else:
        wavemaker_type = args.wavemaker

    if args.mesh == "":
        meshname = "water_no_float0d1"
    else:
        meshname = args.mesh

    id += "_H" + str(H).replace(".", "d")
    id += "_L" + str(L).replace(".", "d")
    id += "_WAVE" + wavemaker_type
    id += "_MESH" + meshname

    if args.dt is not None:
        max_dt = args.dt
    else:
        max_dt = 1 / 20
    id += "_DT" + str(max_dt).replace(".", "d")

    if args.element != "":
        id += "_ELEM" + args.element

    if args.ALE != "":
        id += "_ALE" + args.ALE

    if args.ALEPERIOD != "":
        id += "_ALEPERIOD" + args.ALEPERIOD

    if args.suffix != "":
        id += "_" + args.suffix

    if "piston" in wavemaker_type:
        wavemaker = {"name": "wavemaker", "type": "RigidBody", "velocity": ["piston_using_wave_length", start, a, L, h, 1, 0, 0]}
    elif "flap" in wavemaker_type:
        wavemaker = {"name": "wavemaker", "type": "RigidBody", "velocity": ["flap_using_wave_length", start, a, L, h, h, 0, 1, 0]}
    elif "potential" in wavemaker_type:
        wavemaker = {"name": "wavemaker", "type": "SoftBody", "isFixed": True, "velocity": ["velocity_using_wave_length", start, a, L, h, h]}

    Lx = 0.74
    Lz = 0.415
    d = 0.25
    Ly = 1

    floats = []
    objfolder = os.path.join(os.path.dirname(__file__), "../../../cpp/obj/Tanizawa1996_2d")
    m = 184.3
    Ixx = 1. / 12. * m * (Ly * Ly + Lz * Lz)
    K = 0.266
    Iyy = m * K * K
    Izz = 1. / 12. * m * (Lx * Lx + Ly * Ly)
    float = {"name": "float", "type": "RigidBody", "velocity": "floating", "damping": [1000, 1000, 1000, 1000, 1000, 1000, 0., 1.], "isFixed": [False, True, False, True, False, True], "mass": m}
    float["COM"] = [12., 0., z_surface - d + 0.22]
    float["spring"] = [float["COM"][0], float["COM"][1], float["COM"][2], 51., 1000., 0.]
    float["MOI"] = [1000 * Iyy, Iyy, 1000 * Iyy]
    float["objfile"] = os.path.join(objfolder, "float.obj")
    floats.append(float)

    water["objfile"] = os.path.join(objfolder, meshname + ".obj")

    if "no_float" in id:
        wavemaker["objfile"] = os.path.join(objfolder, "wavemaker_no_float.obj")
        tank["objfile"] = os.path.join(objfolder, "tank_no_float.obj")
    else:
        wavemaker["objfile"] = os.path.join(objfolder, "wavemaker.obj")
        tank["objfile"] = os.path.join(objfolder, "tank.obj")

    absorber["objfile"] = os.path.join(objfolder, "absorber.obj")

    inputfiles = [tank, wavemaker, water, absorber]
    inputfiles += floats

    setting = {"max_dt": max_dt, "end_time_step": 20000, "end_time": 30, "element": args.element, "ALE": args.ALE, "ALEPERIOD": args.ALEPERIOD}

    generate_input_files(inputfiles, setting, lambda id: (f"./input_files/{id}", f"{args.outputdir}/{id}"), id)
