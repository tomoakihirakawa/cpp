import json
import sys
import os

from PySide6.QtWidgets import (QApplication, QWidget, QVBoxLayout, QHBoxLayout, QTableWidget, QTableWidgetItem, QListWidget)
from PySide6.QtCore import Qt

file_dir = './input_files/Cheng2018_H0d04_T1d2_piston_with_float/'

class ExampleApp(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        # Create a horizontal layout to hold the two vertical layouts
        self.mainLayout = QHBoxLayout(self)
        
        self.layoutLeft = QVBoxLayout()
        self.layoutRight = QVBoxLayout()
        self.list_input_json = QListWidget()

        self.tableLeft = QTableWidget()
        self.tableLeft.setRowCount(0)
        self.tableLeft.setColumnCount(2)
        self.tableLeft.setHorizontalHeaderLabels(['Name', 'Value'])
        
        self.tableRight = QTableWidget()
        self.tableRight.setRowCount(0)
        self.tableRight.setColumnCount(2)
        self.tableRight.setHorizontalHeaderLabels(['Name', 'Value'])

        self.layoutLeft.addWidget(self.tableLeft)
        self.layoutRight.addWidget(self.tableRight)
        
        # Add the left and right layouts to the main horizontal layout
        self.mainLayout.addLayout(self.layoutLeft)
        self.mainLayout.addLayout(self.layoutRight)
        self.layoutLeft.addWidget(self.list_input_json)

        # Set the main layout to the widget
        self.setLayout(self.mainLayout)

        # Load the json data
        with open(file_dir + 'setting.json') as f:
            data = json.load(f)
            print(data)
            for key, value in data.items():
                rowCount = self.tableLeft.rowCount()
                self.tableLeft.insertRow(rowCount)
                self.tableLeft.setItem(rowCount, 0, QTableWidgetItem(key))
                self.tableLeft.setItem(rowCount, 1, QTableWidgetItem(str(value)))

                if key == "input_files":
                    for file in value:
                        self.list_input_json.addItem(file)

        # display the selected json on the right table
        self.list_input_json.itemClicked.connect(self.display_selected_json)

    def display_selected_json(self, item):
        self.tableRight.setRowCount(0)
        with open(file_dir + item.text()) as f:
            data = json.load(f)
            for key, value in data.items():
                rowCount = self.tableRight.rowCount()
                self.tableRight.insertRow(rowCount)
                self.tableRight.setItem(rowCount, 0, QTableWidgetItem(key))
                self.tableRight.setItem(rowCount, 1, QTableWidgetItem(str(value)))
                
app = QApplication(sys.argv)
ex = ExampleApp()
ex.show()
sys.exit(app.exec_())
