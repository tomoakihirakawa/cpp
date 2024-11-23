import os
from IOgenerator import generate_input_files, read_obj, Norm3d, Differece3d, Add3d

def generate_input_files(args):
    objfolder = os.path.join(os.path.dirname(__file__), "../../../cpp/obj/Palm2016")

    water = {"name": "water", "type": "Fluid"}
    tank = {"name": "tank", "type": "RigidBody", "isFixed": True}
    wavemaker = {"name": "wavemaker", "type": "SoftBody", "isFixed": True}

    start = 0.
    if args.H is None:
        H = 0.04
    else:
        H = args.H

    a = H / 2
    T = 1.
    h = 0.9
    z_surface = 0.9

    id = args.case
    id += "_H" + str(H).replace(".", "d")
    id += "_T" + str(T).replace(".", "d")

    rho = 1000
    draft = 0.172
    D = 0.515
    vol = math.pi * D ** 2. / 4. * draft
    M = rho * vol
    Ixx = 0.9
    float_bottom = z_surface - draft
    COM = [6, 0, 0.0758 + float_bottom]

    i = 0
    horizontal_length = 1.66
    r = D / 2 + 0.015
    q = i * 2 * pi / 3
    X_fair_leadA = [r * math.cos(q) + COM[0], r * math.sin(q) + COM[1], z_surface]
    X_anchorA = Add3d(X_fair_leadA, [horizontal_length * math.cos(q), horizontal_length * math.sin(q), -z_surface])

    q = (i + 1) * 2 * pi / 3
    X_fair_leadB = [r * math.cos(q) + COM[0], r * math.sin(q) + COM[1], z_surface]
    X_anchorB = Add3d(X_fair_leadB, [horizontal_length * math.cos(q), horizontal_length * math.sin(q), -z_surface])

    q = (i + 2) * 2 * pi / 3
    X_fair_leadC = [r * math.cos(q) + COM[0], r * math.sin(q) + COM[1], z_surface]
    X_anchorC = Add3d(X_fair_leadC, [horizontal_length * math.cos(q), horizontal_length * math.sin(q), -z_surface])

    total_lengthA = Norm3d(Differece3d(X_anchorA, X_fair_leadA))
    total_lengthB = Norm3d(Differece3d(X_anchorB, X_fair_leadB))
    total_lengthC = Norm3d(Differece3d(X_anchorC, X_fair_leadC))

    stiffness = 300 * 10 ** 6
    damp = .5
    density = 0.1447
    n_points = 30
    diam = 4.786 / 1000
    stiffness = stiffness * (math.pi * (diam / 2.) ** 2)

    float = {"name": "float", "type": "RigidBody", "velocity": "floating", "mass": M, "COM": COM, "MOI": [Ixx, Ixx, Ixx * 1000.]}

    id += "_MESH" + args.mesh

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

    if "with_mooring" in id:
        float["mooringA"] = ["mooringA", X_anchorA[0], X_anchorA[1], X_anchorA[2], X_fair_leadA[0], X_fair_leadA[1], X_fair_leadA[2], total_lengthA, n_points, density, stiffness, damp, diam]
        float["mooringB"] = ["mooringB", X_anchorB[0], X_anchorB[1], X_anchorB[2], X_fair_leadB[0], X_fair_leadB[1], X_fair_leadB[2], total_lengthB, n_points, density, stiffness, damp, diam]
        float["mooringC"] = ["mooringC", X_anchorC[0], X_anchorC[1], X_anchorC[2], X_fair_leadC[0], X_fair_leadC[1], X_fair_leadC[2], total_lengthC, n_points, density, stiffness, damp, diam]

    probe1 = {"name": "probe1", "type": "wave gauge", "position": [float["COM"][0] - 0.6, float["COM"][1] + 2., 1.2, float["COM"][0] - 0.6, float["COM"][1] + 2., 0.6]}
    probe2 = {"name": "probe2", "type": "wave gauge", "position": [float["COM"][0], float["COM"][1] + 2., 1.2, float["COM"][0], float["COM"][1] + 2., 0.6]}
    probe3 = {"name": "probe3", "type": "wave gauge", "position": [float["COM"][0] + 0.14, float["COM"][1] + 2., 1.2, float["COM"][0] - 0.14, float["COM"][1] + 2., 0.6]}
    probe4 = {"name": "probe4", "type": "wave gauge", "position": [float["COM"][0] - 0.31, float["COM"][1] + 2., 1.2, float["COM"][0] - 0.31, float["COM"][1] + 2., 0.6]}

    probes = [probe1, probe2, probe3, probe4]

    water["objfile"] = os.path.join(objfolder, args.mesh + ".obj")
    wavemaker["objfile"] = os.path.join(objfolder, "wavemaker_mod.obj")
    tank["objfile"] = os.path.join(objfolder, "tank_mod.obj")
    float["objfile"] = os.path.join(objfolder, "float_mod.obj")

    inputfiles = [tank, wavemaker, water, float]
    inputfiles += probes

    setting = {"max_dt": max_dt, "end_time_step": 100000, "end_time": 30, "element": args.element, "ALE": args.ALE, "ALEPERIOD": args.ALEPERIOD}

    generate_input_files(inputfiles, setting, lambda id: (f"./input_files/{id}", f"{args.outputdir}/{id}"), id)
