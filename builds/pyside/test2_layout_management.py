import sys
from PySide6.QtWidgets import QApplication, QMainWindow, QWidget, QVBoxLayout, QPushButton

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Layout Management")
        layout = QVBoxLayout()
        button1 = QPushButton("Button 1")
        layout.addWidget(button1)
        button1.clicked.connect(self.on_button_clicked)
        self.button1_clicked = False
        layout.addWidget(QPushButton("Button 2"))
        layout.addWidget(QPushButton("Button 3"))

        container = QWidget()
        container.setLayout(layout)

        self.setCentralWidget(container)
        
    def on_button_clicked(self):
        self.button1_clicked = not self.button1_clicked
        # change button color
        if self.button1_clicked:
            self.sender().setStyleSheet("background-color: lightblue")
        else:
            self.sender().setStyleSheet("background-color: lightgreen")


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())