def generate_input_files(args):
    start = 0.

    id = args.SimulationCase

    objfolder = args.code_home_dir + "/cpp/obj/fish"

    waterA = {"name": "waterA",
              "type": "Fluid",
              "reverseNormal": True,
              "objfile": objfolder + "/bodyA20.obj"}

    waterB = {"name": "waterB",
              "type": "Fluid",
              "reverseNormal": True,
              "objfile": objfolder + "/bodyB20.obj"}

    waterC = {"name": "waterC",
              "type": "Fluid",
              "reverseNormal": True,
              "objfile": objfolder + "/bodyC20.obj"}

    bodyA = {"name": "bodyA",
             "type": "RigidBody",
             "COM": [0., 0, 0.25],
             "mass": 10**10,
             "MOI": [10**10, 10**10, 10**10],
             "output": "json",
             "velocity": ["file", "./study_fish/bodyA.dat"],
             "objfile": objfolder + "/bodyA20.obj"}

    bodyB = {"name": "bodyB",
             "type": "RigidBody",
             "COM": [0.35, 0, 0.25],
             "mass": 10**10,
             "MOI": [10**10, 10**10, 10**10],
             "output": "json",
             "velocity": ["file", "./study_fish/bodyB.dat"],
             "objfile": objfolder + "/bodyB20.obj"}

    bodyC = {"name": "bodyC",
             "type": "RigidBody",
             "COM": [0.7, 0, 0.25],
             "mass": 10**10,
             "MOI": [10**10, 10**10, 10**10],
             "output": "json",
             "velocity": ["file", "./study_fish/bodyC.dat"],
             "objfile": objfolder + "/bodyC20.obj"}

    inputfiles = [waterA, waterB, waterC, bodyA, bodyB, bodyC]

    setting = {"max_dt": 0.01,
               "end_time_step": 10000,
               "end_time": 9}

    args.generate_input_files(inputfiles, setting, args.IO_dir, id)
