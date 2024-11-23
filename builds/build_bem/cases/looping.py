import os
from IOgenerator import generate_input_files, read_obj

def generate_input_files(args):
    objfolder = os.path.join(os.path.dirname(__file__), "../../../cpp/obj/looping")

    water = {"name": "water", "type": "Fluid", "objfile": os.path.join(objfolder, "water_remesh.obj")}
    vertices, triangles = read_obj(water["objfile"])
    water["vertices"] = vertices
    water["triangles"] = triangles

    tank = {"name": "tank", "type": "RigidBody", "isFixed": True, "objfile": os.path.join(objfolder, "tank.obj")}

    start = 0.
    T = 1.
    H = 0.05
    a = H / 2
    h = 1.5
    z_surface = h

    id = args.case
    id += "_H" + str(H).replace(".", "d")
    id += "_T" + str(T).replace(".", "d")
    id += "_h" + str(h).replace(".", "d")

    if "piston" in args.case:
        wavemaker = {"name": "wavemaker", "type": "RigidBody", "velocity": ["piston", start, a, T, h, 1, 0, 0], "objfile": os.path.join(objfolder, "wavemaker.obj")}
    elif "flap" in args.case:
        wavemaker = {"name": "wavemaker", "type": "RigidBody", "velocity": ["flap", start, a, T, h, h, 0, 1, 0], "objfile": os.path.join(objfolder, "wavemaker.obj")}
    else:
        wavemaker = {"name": "wavemaker", "type": "SoftBody", "isFixed": True, "velocity": ["velocity", start, a, T, h, z_surface], "objfile": os.path.join(objfolder, "wavemaker.obj")}

    inputfiles = [tank, wavemaker, water]

    setting = {"max_dt": 0.03, "end_time_step": 10000, "end_time": 9}

    generate_input_files(inputfiles, setting, lambda id: (f"./input_files/{id}", f"{args.outputdir}/{id}"), id)
