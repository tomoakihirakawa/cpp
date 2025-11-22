# check_paths.py
import os

# このスクリプト(check_paths.py)があるディレクトリを取得
current_dir = os.path.dirname(os.path.abspath(__file__))
print(f"スクリプトの場所: {current_dir}")

# 3階層上のディレクトリの絶対パスを計算して表示
# (あなたのフォルダ構成に合わせて ../ の数を調整してください)
project_root = os.path.abspath(os.path.join(current_dir, "..", ".."))
print(f"\n想定されるプロジェクトルート: {project_root}")

# 想定されるobjフォルダのパス
obj_path = os.path.join(project_root, "cpp", "obj")
print(f"\nこのパスをコピーしてください↓\n{obj_path}")