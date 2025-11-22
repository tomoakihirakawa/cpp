# input_generator_logic.py
import os

# ▼▼▼ ここに、check_paths.pyの出力結果（最後の行）を貼り付け ▼▼▼
# 例：BASE_OBJ_PATH = "/Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/obj"
BASE_OBJ_PATH = "/Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/obj"
# ▲▲▲ 必ず書き換えてください ▲▲▲


def get_case_data(case_name, args_dict={}, base_dir="."):
    """
    指定されたケース名と引数に基づいて、シミュレーションの設定データ（辞書）を生成して返す。
    """
    # --- この関数の中身は、BASE_OBJ_PATH を使うように変更 ---

    dt = args_dict.get('max_dt', 0.05)
    suffix = args_dict.get('suffix', 'Trial10')

    if "Ruehl2016" in case_name:
        objfolder = os.path.join(BASE_OBJ_PATH, "Ruehl2016")

        # study_benchmarkフォルダへのパスも修正
        # BASE_OBJ_PATHから逆算してプロジェクトルートを特定
        project_root = os.path.abspath(os.path.join(BASE_OBJ_PATH, "..", ".."))
        benchmark_folder = os.path.join(project_root, "study_benchmark", "Ruehl2016")

        h = 1.36
        trial = suffix.split("_")[0] if "_" in suffix else suffix

        water = {"name": "water", "type": "Fluid", "objfile": os.path.join(objfolder, "water.obj")}
        tank = {"name": "tank", "type": "RigidBody", "isFixed": True, "objfile": os.path.join(objfolder, "tank.obj")}
        wavemaker = {"name": "wavemaker", "type": "RigidBody", "velocity": ["file", os.path.join(benchmark_folder, f"{trial}.csv")], "objfile": os.path.join(objfolder, "wavemaker.obj")}
        absorber = {"name": "absorber", "type": "Absorber", "isFixed": True, "objfile": os.path.join(objfolder, "absorber.obj"), "wave_theory_L": [0, 1., h, 0]}
        wg3 = {"name": "wg3", "type": "wave gauge", "position": [9.48, 0., 1.36+0.5, 9.48, 0., 1.36-0.5]}
        wg6 = {"name": "wg6", "type": "wave gauge", "position": [16.778, -0.031, 2., 16.778, -0.031, 0.5]}

        inputfiles = [tank, wavemaker, water, absorber, wg3, wg6]
        setting = { "max_dt": dt, "end_time": 35, "element": "linear", "ALE": "linear", "ALEPERIOD": "1" }

        return setting, inputfiles

    elif "Tanizawa1996" in case_name:
        objfolder = os.path.join(BASE_OBJ_PATH, "Tanizawa1996")

        water = {"name": "water", "type": "Fluid", "objfile": os.path.join(objfolder, "water.obj")}
        tank = {"name": "tank", "type": "RigidBody", "isFixed": True, "objfile": os.path.join(objfolder, "tank.obj")}
        box = {"name": "floating_box", "type": "RigidBody", "isFixed": False, "objfile": os.path.join(objfolder, "box.obj"), "mass": 100, "COM": [0,0,0]}

        inputfiles = [water, tank, box]
        setting = { "max_dt": 0.01, "end_time": 20.0 }

        return setting, inputfiles

    return None, None