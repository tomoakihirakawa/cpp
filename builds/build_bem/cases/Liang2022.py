import os
from IOgenerator import generate_input_files, read_obj, Norm3d, Differece3d, Add3d

def generate_input_files(args):
    objfolder = os.path.join(os.path.dirname(__file__), "../../../cpp/obj/Liang2022")

    water = {"name": "water", "type": "Fluid"}
    tank = {"name": "tank", "type": "RigidBody", "isFixed": True}
    wavemaker = {"name": "wavemaker", "type": "RigidBody"}

    start = 0.
    if args.H is None:
        H = 0.07
    else:
        H = args.H

    a = H / 2
    T = 1.4
    h = 0.6
    z_surface = h

    id = args.case
    id += "_H" + str(H).replace(".", "d")
    id += "_T" + str(T).replace(".", "d")

    rho = 1000
    draft = 0.16
    Lx = W = 0.5  # float width in x direction
    Ly = 0.745

    V = draft * Lx * Ly
    print("V ", V, ", mass should be ", V * rho)

    M = 58.09  # kg (質量)
    Iyy = 2.441  # kg*m^2 (慣性モーメント)
    float_bottom = z_surface - draft  # m (底面の高さ)
    COM = [4.5 + 0.5 / 2, 0, 0.0652 + float_bottom]  # m (重心位置)

    float = {"name": "float", "type": "RigidBody", "velocity": "floating", "mass": M, "COM": COM, "MOI": [10 ** 10, Iyy, 10 ** 10]}

    if "typeA" in args.case:
        r = math.sqrt(1.567 ** 2 - float_bottom ** 2)  # horizontal_length
        X_fair_leadA = [float["COM"][0] - Lx / 2, float["COM"][1] - Ly / 2, float_bottom]
        X_anchorA = [float["COM"][0] - Lx / 2 + r, float["COM"][1] - Ly / 2, 0]
        X_fair_leadB = [float["COM"][0] - Lx / 2, float["COM"][1] + Ly / 2, float_bottom]
        X_anchorB = [float["COM"][0] - Lx / 2 + r, float["COM"][1] + Ly / 2, 0]
        X_fair_leadC = [float["COM"][0] + Lx / 2, float["COM"][1] - Ly / 2, float_bottom]
        X_anchorC = [float["COM"][0] + Lx / 2 - r, float["COM"][1] - Ly / 2, 0]
        X_fair_leadD = [float["COM"][0] + Lx / 2, float["COM"][1] + Ly / 2, float_bottom]
        X_anchorD = [float["COM"][0] + Lx / 2 - r, float["COM"][1] + Ly / 2, 0]
    elif "typeB" in args.case:
        r = math.sqrt(1.125 ** 2 - float_bottom ** 2)
        X_fair_leadA = [float["COM"][0] - Lx / 2, float["COM"][1] - Ly / 2, float_bottom]
        X_anchorA = [float["COM"][0] - Lx / 2 - r, float["COM"][1] - Ly / 2, 0]
        X_fair_leadB = [float["COM"][0] - Lx / 2, float["COM"][1] + Ly / 2, float_bottom]
        X_anchorB = [float["COM"][0] - Lx / 2 - r, float["COM"][1] + Ly / 2, 0]
        X_fair_leadC = [float["COM"][0] + Lx / 2, float["COM"][1] - Ly / 2, float_bottom]
        X_anchorC = [float["COM"][0] + Lx / 2 + r, float["COM"][1] - Ly / 2, 0]
        X_fair_leadD = [float["COM"][0] + Lx / 2, float["COM"][1] + Ly / 2, float_bottom]
        X_anchorD = [float["COM"][0] + Lx / 2 + r, float["COM"][1] + Ly / 2, 0]
    elif "typeC" in args.case:
        r = math.sqrt(0.809 ** 2 - float_bottom ** 2)
        X_fair_leadA = [float["COM"][0] - Lx / 2, float["COM"][1] - Ly / 2, float_bottom]
        X_anchorA = [float["COM"][0] - Lx / 2 - r, float["COM"][1] - Ly / 2, 0]
        X_fair_leadB = [float["COM"][0] - Lx / 2, float["COM"][1] + Ly / 2, float_bottom]
        X_anchorB = [float["COM"][0] - Lx / 2 - r, float["COM"][1] + Ly / 2, 0]
        X_fair_leadC = [float["COM"][0] + Lx / 2, float["COM"][1] - Ly / 2, float_bottom]
        X_anchorC = [float["COM"][0] + Lx / 2 + r, float["COM"][1] - Ly / 2, 0]
        X_fair_leadD = [float["COM"][0] + Lx / 2, float["COM"][1] + Ly / 2, float_bottom]
        X_anchorD = [float["COM"][0] + Lx / 2 + r, float["COM"][1] + Ly / 2, 0]

    stiffness = 2.36 * 10 ** 3  # ! [N/m]
    damp = .4  # ! [N/(m/s^2)]
    density = 0.177  # ! [kg/m]

    total_lengthA = Norm3d(Differece3d(X_anchorA, X_fair_leadA))
    total_lengthB = Norm3d(Differece3d(X_anchorB, X_fair_leadB))
    total_lengthC = Norm3d(Differece3d(X_anchorC, X_fair_leadC))
    total_lengthD = Norm3d(Differece3d(X_anchorD, X_fair_leadD))

    n_points = 30
    diam = 1 * 0.1  # 1 cm

    float["mooringA"] = ["mooringA", X_anchorA[0], X_anchorA[1], X_anchorA[2], X_fair_leadA[0], X_fair_leadA[1], X_fair_leadA[2], total_lengthA, n_points, density, stiffness, damp, diam]
    float["mooringB"] = ["mooringB", X_anchorB[0], X_anchorB[1], X_anchorB[2], X_fair_leadB[0], X_fair_leadB[1], X_fair_leadB[2], total_lengthB, n_points, density, stiffness, damp, diam]
    float["mooringC"] = ["mooringC", X_anchorC[0], X_anchorC[1], X_anchorC[2], X_fair_leadC[0], X_fair_leadC[1], X_fair_leadC[2], total_lengthC, n_points, density, stiffness, damp, diam]
    float["mooringD"] = ["mooringD", X_anchorD[0], X_anchorD[1], X_anchorD[2], X_fair_leadD[0], X_fair_leadD[1], X_fair_leadD[2], total_lengthD, n_points, density, stiffness, damp, diam]

    probe1 = {"name": "WG1", "type": "wave gauge", "position": [float["COM"][0] - W / 2 - 3.5, float["COM"][1], 1.2, float["COM"][0] - W / 2 - 3.5, float["COM"][1], 0.3]}
    probe2 = {"name": "WG2", "type": "wave gauge", "position": [float["COM"][0] - W / 2 - 3.0, float["COM"][1], 1.2, float["COM"][0] - W / 2 - 3.0, float["COM"][1], 0.3]}
    probe3 = {"name": "WG3", "type": "wave gauge", "position": [float["COM"][0] + W / 2 + 3.0, float["COM"][1], 1.2, float["COM"][0] + W / 2 + 3.0, float["COM"][1], 0.3]}
    probe4 = {"name": "WG4", "type": "wave gauge", "position": [float["COM"][0] + W / 2 + 3.5, float["COM"][1], 1.2, float["COM"][0] + W / 2 + 3.5, float["COM"][1], 0.3]}
    probes = [probe1, probe2, probe3, probe4]

    water["objfile"] = os.path.join(objfolder, "water5_mod.obj")
    wavemaker["objfile"] = os.path.join(objfolder, "wavemaker5.obj")
    tank["objfile"] = os.path.join(objfolder, "tank5.obj")
    float["objfile"] = os.path.join(objfolder, "float5.obj")

    inputfiles = [tank, wavemaker, water, float]
    inputfiles += probes

    setting = {"max_dt": 0.03, "end_time_step": 100000, "end_time": 30}

    generate_input_files(inputfiles, setting, lambda id: (f"./input_files/{id}", f"{args.outputdir}/{id}"), id)
