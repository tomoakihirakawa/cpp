
import json 
import sys
from PySide6.QtWidgets import (QApplication, QWidget, QVBoxLayout, QTableWidget, QTableWidgetItem, QPushButton, QComboBox, QLineEdit, QLabel, QMenu)
from PySide6.QtCore import Qt

class ExampleApp(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        self.layout = QVBoxLayout(self)
        self.tableWidget = QTableWidget()
        self.tableWidget.setRowCount(0)
        self.tableWidget.setColumnCount(2)
        self.tableWidget.setHorizontalHeaderLabels(['Name', 'Value'])
        self.layout.addWidget(self.tableWidget)
        self.setLayout(self.layout)

        # load the json data
        with open('./input_files/Cheng2018_free_decay_H0d04_T1d2_piston_with_float/setting.json') as f:
            data = json.load(f)
            print(data)
            for key, value in data.items():
                self.tableWidget.insertRow(self.tableWidget.rowCount())
                self.tableWidget.setItem(self.tableWidget.rowCount() - 1, 0, QTableWidgetItem(key))
                self.tableWidget.setItem(self.tableWidget.rowCount() - 1, 1, QTableWidgetItem(str(value)))


app = QApplication(sys.argv)
ex = ExampleApp()
ex.show()
sys.exit(app.exec_())
