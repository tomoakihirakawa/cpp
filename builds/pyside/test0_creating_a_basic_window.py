import sys
from PySide6.QtWidgets import QApplication, QMainWindow

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("PySide6 Basic Window")
        self.setGeometry(100,100,400,300)
        # display the window size with arrow
        self.statusBar().showMessage("Size: 400x300")
        self.statusBar().setStyleSheet("background-color: lightblue")
        self.statusBar().setStyleSheet("color: red")
        


def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()