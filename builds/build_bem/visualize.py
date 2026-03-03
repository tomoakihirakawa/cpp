import sys
import json
import os

try:
    from paraview.simple import *
except ImportError:
    print("エラー: ParaViewライブラリが見つかりません。pvpythonで実行してください。")
    sys.exit(1)

# ★★★ ここからが変更点 ★★★
# オブジェクト名と色(RGB)を対応させる辞書を定義
COLOR_PALETTE = {
    "water":     [0.3, 0.5, 0.9],  # 青
    "tank":      [0.6, 0.6, 0.6],  # グレー
    "float":     [1.0, 0.84, 0.0], # 黄色
    "wavemaker": [0.2, 0.8, 0.2],  # 緑
    "absorber":  [0.8, 0.5, 0.2],  # オレンジ
}
# ★★★ ここまでが変更点 ★★★

# --- 環境変数からIDを取得 ---
try:
    simulation_id_or_path = os.environ.get('SIM_ID')
    if not simulation_id_or_path:
        raise ValueError("エラー: 環境変数 'SIM_ID' が設定されていません。launch_paraview.sh を使って実行してください。")
except Exception as e:
    try:
        from paraview.simple import Show, Text
        Show(Text(Text=str(e)))
    except:
        pass
    raise

# --- 入力ディレクトリのパスを決定 ---
if "input_files" in simulation_id_or_path:
    input_dir = simulation_id_or_path
else:
    input_dir = os.path.join("./input_files", simulation_id_or_path)
print(f"入力ディレクトリ '{input_dir}' のシーンを構築します...")

# --- settings.json を読み込む ---
try:
    settings_path = os.path.join(input_dir, "settings.json")
    if not os.path.exists(settings_path):
        # Backward compatibility
        settings_path = os.path.join(input_dir, "setting.json")
    with open(settings_path) as f:
        settings = json.load(f)
except FileNotFoundError:
    print(f"エラー: {input_dir} に settings.json が見つかりません。")
    raise

# --- ParaViewのシーンを構築 ---
renderView = GetActiveViewOrCreate('RenderView')
renderView.ViewSize = [1200, 800]
renderView.Background = [0.1, 0.1, 0.15] # 背景を濃いグレーに

# --- メインの処理ループ ---
for file_name in settings.get("input_files", []):
    try:
        with open(os.path.join(input_dir, file_name)) as f:
            info = json.load(f)
            object_name = info.get("name", "Unnamed").lower() # 小文字に統一
            object_type = info.get("type", "").lower()

            # ケース1: .objファイルを持つオブジェクト
            if "objfile" in info and os.path.exists(info["objfile"]):
                reader = OpenDataFile(info["objfile"])
                display = Show(reader, renderView)
                display.Representation = 'Surface'
                RenameSource(info.get("name"), reader)

                # ★★★ ここからが変更点 ★★★
                # 名前に応じて色と透明度を設定
                if "water" in object_name:
                    display.DiffuseColor = COLOR_PALETTE["water"]
                    display.Opacity = 0.4 # 水を半透明に
                elif "tank" in object_name or "wall" in object_name:
                    display.DiffuseColor = COLOR_PALETTE["tank"]
                elif "float" in object_name:
                    display.DiffuseColor = COLOR_PALETTE["float"]
                elif "wavemaker" in object_name or "maker" in object_name:
                    display.DiffuseColor = COLOR_PALETTE["wavemaker"]
                elif "absorber" in object_name:
                    display.DiffuseColor = COLOR_PALETTE["absorber"]
                # ★★★ ここまでが変更点 ★★★

            # ケース2: 波高計
            elif "wave gauge" in object_type:
                if "position" in info and len(info["position"]) == 6:
                    pos = info["position"]
                    p1, p2 = [pos[0], pos[1], pos[2]], [pos[3], pos[4], pos[5]]
                    line = LineSource(Point1=p1, Point2=p2)
                    display = Show(line, renderView)
                    display.DiffuseColor = [1.0, 0.2, 0.2] # 赤色
                    display.LineWidth = 4.0
                    RenameSource(info.get("name"), line)

    except (FileNotFoundError, json.JSONDecodeError) as e:
        print(f"警告: {file_name} の処理中に問題が発生しました - {e}")
        pass
    except Exception as e:
        print(f"処理中に予期せぬエラーが発生しました: {e}")


# --- 設定情報をテキストで表示 ---
text_source = Text(Text=json.dumps(settings, indent=2, ensure_ascii=False))
text_display = Show(text_source, renderView)
text_display.WindowLocation = 'Lower Left Corner'
text_display.FontSize = 8
RenameSource("Settings", text_source)

# --- カメラをリセットして全体を表示 ---
Render()
ResetCamera()
