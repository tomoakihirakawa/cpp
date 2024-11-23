def generate_input_files(args):
    start = 0.

    id = args.SimulationCase

    objfolder = args.code_home_dir + "/cpp/obj/fish"
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
             "velocity": ["file", "./study_fish/bodyA.dat"],
             "objfile": objfolder + "/bodyA50.obj"}

    bodyB = {"name": "bodyB",
             "type": "RigidBody",
             "COM": [0.35, 0, 0.25],
             "mass": 10**10,
             "MOI": [10**10, 10**10, 10**10],
             "output": "json",
             "velocity": ["file", "./study_fish/bodyB.dat"],
             "objfile": objfolder + "/bodyB50.obj"}

    bodyC = {"name": "bodyC",
             "type": "RigidBody",
             "COM": [0.7, 0, 0.25],
             "mass": 10**10,
             "MOI": [10**10, 10**10, 10**10],
             "output": "json",
             "velocity": ["file", "./study_fish/bodyC.dat"],
             "objfile": objfolder + "/bodyC50.obj"}

    inputfiles = [tank, water, bodyA, bodyB, bodyC]

    setting = {"max_dt": 0.01,
               "end_time_step": 10000,
               "end_time": 9}

    args.generate_input_files(inputfiles, setting, args.IO_dir, id)
