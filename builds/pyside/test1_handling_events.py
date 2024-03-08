'''DOC_EXTRACT

`QPushButton`には，`clicked`というシグナルがあります．
`clicked`には，ボタンがクリックされたときに実行される関数を接続`connect`することができる．

```python
button.clicked.connect(self.on_button_clicked)
```

'''

import sys
from PySide6.QtWidgets import QApplication, QMainWindow, QPushButton

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Event Handling")
        button = QPushButton("Click me for a message!")
        button.clicked.connect(self.on_button_clicked)
        self.setCentralWidget(button)
        self.switch = False
        

    def on_button_clicked(self):
        # change color of the button
        self.switch = not self.switch
        if self.switch:
            self.sender().setStyleSheet("background-color: lightblue")
        else:
            self.sender().setStyleSheet("background-color: lightgreen") 


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
