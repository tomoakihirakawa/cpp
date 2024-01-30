import sys
from PySide6.QtWidgets import QApplication, QWidget, QVBoxLayout, QLineEdit

app = QApplication(sys.argv)
window = QWidget()
layout = QVBoxLayout()

edit = QLineEdit()
layout.addWidget(edit)

window.setLayout(layout)
window.show()
sys.exit(app.exec_())
