import os
from IOgenerator import generate_input_files, read_obj

def generate_input_files(args):
    objfolder = os.path.join(os.path.dirname(__file__), "../../../cpp/obj/Hadzic2005")

    water = {"name": "water", "type": "Fluid"}
    if "multiple" in args.case:
        water["objfile"] = os.path.join(objfolder, "water1000meshlab.obj")
    elif "without_float" in args.case:
        water["objfile"] = os.path.join(objfolder, "water_without_float9.obj")
    else:
        water["objfile"] = os.path.join(objfolder, "water1000.obj")

    vertices, triangles = read_obj(water["objfile"])
    water["vertices"] = vertices
    water["triangles"] = triangles

    tank = {"name": "tank", "type": "RigidBody", "isFixed": True, "objfile": os.path.join(objfolder, "tank10.obj")}

    start_time = 0.
    wavemaker = {"name": "wavemaker", "type": "RigidBody", "objfile": os.path.join(objfolder, "wavemaker100.obj"), "velocity": ["Hadzic2005", start_time], "acceleration": ["Hadzic2005", start_time], "COM": [0., 0., 0.]}

    float = {"name": "float", "type": "RigidBody", "objfile": os.path.join(objfolder, "float20.obj"), "output": "json", "velocity": "floating"}

    L = 0.1
    W = 0.29
    H = 0.05
    A = L * W
    density = 680
    float["mass"] = m = density * 0.05 * 0.1 * 0.29
    d = m / (1000 * A)
    Ixx = 1. / 12. * m * (W * W + H * H)
    Iyy = 1. / 12. * m * (L * L + H * H)
    Izz = 1. / 12. * m * (L * L + W * W)
    z_surface = 0.4
    z_floatinbody_bottom = z_surface - d
    float["COM"] = [2.11, W / 2, z_floatinbody_bottom + H / 2]
    float["MOI"] = [Ixx * 10 ** 10, 0.001589233395, Izz * 10 ** 10]

    if "without_float" in args.case:
        inputfiles = [tank, wavemaker, water]
    else:
        inputfiles = [tank, wavemaker, water, float]

    setting = {"max_dt": 0.05, "end_time_step": 10000, "end_time": 9}

    generate_input_files(inputfiles, setting, lambda id: (f"./input_files/{id}", f"{args.outputdir}/{id}"), args.case)