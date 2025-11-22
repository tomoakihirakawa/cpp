import sys
from PySide6.QtWidgets import QApplication, QMainWindow, QPushButton

def on_button_clicked():
    print("ボタンがクリックされました！")

app = QApplication(sys.argv)
window = QMainWindow()
button = QPushButton("クリックしてください", window)
button.clicked.connect(on_button_clicked)
window.show()
sys.exit(app.exec_())
