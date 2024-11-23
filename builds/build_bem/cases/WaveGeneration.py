def generate_input_files(args):
    objfolder = args.code_home_dir + "/cpp/obj/WaveGeneration"

    water = {"name": "water",
             "type": "Fluid",
             "objfile": objfolder + "/water_without_float9_mod.obj"}

    vertices, triangles = args.read_obj(water["objfile"])
    water["vertices"] = vertices
    water["triangles"] = triangles

    tank = {"name": "tank", "type":
            "RigidBody", "isFixed": True,
            "objfile": objfolder + "/tank10.obj"}

    start = 0.

    T = 1.
    H = 0.05
    a = H/2
    h = 0.4

    id = args.SimulationCase
    id += "_H" + str(H).replace(".", "d")
    id += "_T"+str(T).replace(".", "d")
    id += "_h"+str(h).replace(".", "d")

    z_surface = h

    if "piston" in args.SimulationCase:
        wavemaker = {"name": "wavemaker",
                     "type": "RigidBody",
                     "velocity": ["piston", start, a, T, h, 1, 0, 0],
                     "objfile": f"{objfolder}/wavemaker10.obj"}

    elif "flap" in args.SimulationCase:
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

    args.generate_input_files(inputfiles, setting, args.IO_dir, id)
