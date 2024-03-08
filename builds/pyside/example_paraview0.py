import subprocess

def run_paraview_script():
    # ParaView スクリプトのパス
    script_path = "script.py"
    
    # pvpython を使用してスクリプトを実行
    subprocess.run(["pvpython", script_path])

from PySide6.QtWidgets import QApplication, QPushButton, QVBoxLayout, QWidget

def run_paraview_script():
    # ここにParaViewスクリプトを実行するコードを追加します。
    # 例えば、subprocessを使用して外部のpvpythonスクリプトを実行することができます。
    print("ParaViewスクリプトを実行")

app = QApplication([])

# GUIのレイアウトを設定
window = QWidget()
layout = QVBoxLayout()

# ボタンの追加
button = QPushButton("ParaView スクリプトを実行")
button.clicked.connect(run_paraview_script)
layout.addWidget(button)

window.setLayout(layout)
window.show()

app.exec()
