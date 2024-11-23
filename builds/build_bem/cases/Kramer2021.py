import os
from IOgenerator import generate_input_files, read_obj

def generate_input_files(args):
    objfolder = os.path.join(os.path.dirname(__file__), "../../../cpp/obj/Kramer2021")

    water = {"name": "water", "type": "Fluid"}
    tank = {"name": "tank", "type": "RigidBody", "isFixed": True}

    start_t = 0.02
    float = {"name": "float", "type": "RigidBody", "velocity": ["floating", start_t]}

    float["mass"] = m = 7.056
    z_surface = 900 / 1000
    float["COM"] = [0., 0., z_surface + 0.03]  # 今回は重要ではない
    float["radius_of_gyration"] = [10 ** 10, 10 ** 10, 10 ** 10]
    float["MOI"] = [m * math.pow(float["radius_of_gyration"][0], 2), m * math.pow(float["radius_of_gyration"][1], 2), m * math.pow(float["radius_of_gyration"][2], 2)]

    mesh = "meshB"
    id = args.case + "_H00d03_" + mesh

    if "H00d03" in id:
        objfolder = os.path.join(os.path.dirname(__file__), "../../../cpp/obj/Kramer2021_01D")
        water["objfile"] = os.path.join(objfolder, f"water5_{mesh}.obj")
    elif "H00d09" in id:
        objfolder = os.path.join(os.path.dirname(__file__), "../../../cpp/obj/Kramer2021_03D")
        water["objfile"] = os.path.join(objfolder, f"water4_{mesh}.obj")

    tank["objfile"] = os.path.join(objfolder, "tank4.obj")
    float["objfile"] = os.path.join(objfolder, "sphere.obj")

    gauges = []
    dx = 0.6
    for i in range(10):
        gauges.append({"name": f"gauge{i}", "type": "wave gauge", "position": [0. + dx * i, 0, 1.1, 0. + dx * i, 0, 1.1 - 0.6]})

    inputfiles = [tank, water, float] + gauges

    setting = {
        "WATER_DENSITY": 998.2,
        "GRAVITY": 9.82,
        "max_dt": 0.001,
        "end_time_step": 10000,
        "end_time": 4,
        "ALE_period_in_step": 1
    }

    generate_input_files(inputfiles, setting, lambda id: (f"./input_files/{id}", f"{args.outputdir}/{id}"), id)
