import json
import sys
from PySide6.QtWidgets import (QApplication, QWidget, QVBoxLayout, QLineEdit, QPushButton, 
                               QTableWidget, QTableWidgetItem, QComboBox, QHBoxLayout, QLabel, QMenu)

from PySide6.QtCore import Qt

class ExampleApp(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        self.layout = QVBoxLayout(self)

        # 入力フィールドとコンボボックスのレイアウト
        input_layout = QHBoxLayout()

        # COM ボタンと3Dベクトルの値ボックスの初期化（非表示）
        self.comButton = QPushButton("COM", self)
        self.comButton.setVisible(False)
        self.layout.addWidget(self.comButton)

        self.vectorLayout = QHBoxLayout()
        self.vectorX = QLineEdit(self)
        self.vectorY = QLineEdit(self)
        self.vectorZ = QLineEdit(self)
        self.vectorLayout.addWidget(QLabel("X:"))
        self.vectorLayout.addWidget(self.vectorX)
        self.vectorLayout.addWidget(QLabel("Y:"))
        self.vectorLayout.addWidget(self.vectorY)
        self.vectorLayout.addWidget(QLabel("Z:"))
        self.vectorLayout.addWidget(self.vectorZ)
        self.vectorX.setVisible(False)
        self.vectorY.setVisible(False)
        self.vectorZ.setVisible(False)
        self.layout.addLayout(self.vectorLayout)

        # 名前入力フィールド
        self.nameEdit = QLineEdit(self)
        input_layout.addWidget(QLabel("Name:"))
        input_layout.addWidget(self.nameEdit)

        # タイプ選択コンボボックス
        self.typeCombo = QComboBox(self)
        self.typeCombo.addItems(["water", "rigid body", "softbody"])
        input_layout.addWidget(QLabel("Type:"))
        input_layout.addWidget(self.typeCombo)
        
        # Connect the signal after defining self.typeCombo
        self.typeCombo.currentIndexChanged.connect(self.typeChanged)

        self.layout.addLayout(input_layout)

        # 追加ボタン
        self.addButton = QPushButton("追加", self)
        self.addButton.clicked.connect(self.addEntryToTable)
        self.layout.addWidget(self.addButton)

        # 出力ボタン
        self.outputButton = QPushButton("出力", self)
        self.outputButton.clicked.connect(self.outputData)
        self.layout.addWidget(self.outputButton)

        # テーブルウィジェット
        self.tableWidget = QTableWidget(self)
        self.tableWidget.setColumnCount(2)
        self.tableWidget.setHorizontalHeaderLabels(["Name", "Type"])
        self.layout.addWidget(self.tableWidget)

        self.setLayout(self.layout)

        # コンテキストメニューの設定
        self.tableWidget.setContextMenuPolicy(Qt.CustomContextMenu)
        self.tableWidget.customContextMenuRequested.connect(self.openMenu)

        self.layout.addWidget(self.tableWidget)


    def typeChanged(self, index):
        # "rigid body"が選択された場合にのみCOMボタンとベクトルフィールドを表示
        if self.typeCombo.currentText() == "rigid body":
            self.comButton.setVisible(True)
            self.vectorX.setVisible(True)
            self.vectorY.setVisible(True)
            self.vectorZ.setVisible(True)
        else:
            self.comButton.setVisible(False)
            self.vectorX.setVisible(False)
            self.vectorY.setVisible(False)
            self.vectorZ.setVisible(False)

    def addEntryToTable(self):
        name = self.nameEdit.text()
        type = self.typeCombo.currentText()

        # 名前とタイプが両方入力されているか確認
        if name and type:
            row_count = self.tableWidget.rowCount()
            self.tableWidget.insertRow(row_count)
            self.tableWidget.setItem(row_count, 0, QTableWidgetItem(name))
            self.tableWidget.setItem(row_count, 1, QTableWidgetItem(type))

            # フィールドをクリア
            self.nameEdit.clear()
        else:
            # エラーメッセージを表示するか、入力を促す
            print("名前とタイプを入力してください。")



    def openMenu(self, position):
        menu = QMenu()

        deleteAction = menu.addAction("削除")
        action = menu.exec_(self.tableWidget.mapToGlobal(position))

        if action == deleteAction:
            self.deleteRow()

    def deleteRow(self):
        current_row = self.tableWidget.currentRow()
        if current_row >= 0:
            self.tableWidget.removeRow(current_row)

    def outputData(self):
        for row in range(self.tableWidget.rowCount()):
            name = self.tableWidget.item(row, 0).text()
            type = self.tableWidget.item(row, 1).text()

            data = {'name': name, 'type': type}
            with open(f'{name}.json', 'w') as outfile:
                json.dump(data, outfile)

app = QApplication(sys.argv)
ex = ExampleApp()
ex.show()
sys.exit(app.exec_())