import os
import re
from collections import defaultdict
import datetime 
import os
import json
import sys
from PySide6.QtWidgets import QApplication, QWidget, QLineEdit, QFormLayout, QTreeView, QVBoxLayout, QMessageBox, QHeaderView, QSizePolicy, QFileDialog, QHBoxLayout, QPushButton, QListWidget, QFileSystemModel, QSplitter, QLabel, QTabWidget, QTextEdit
from TreeJSONUtilities import JsonModel  # Assuming TreeItem is used internally in JsonModel
from PySide6.QtCore import Qt, QDir
from PySide6.QtGui import QMouseEvent, QPalette, QColor, QFont

class ClickableLineEdit(QLineEdit):
    def __init__(self, parent=None):
        super().__init__(parent)
    
    def mousePressEvent(self, event: QMouseEvent):
        super().mousePressEvent(event)
        if event.button() == Qt.LeftButton:
            self.openFileDialog()
    
    def openFileDialog(self):
        dirPath = QFileDialog.getExistingDirectory(self, "Select Directory")
        if dirPath:
            self.setText(dirPath)

class MainWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        self.setWindowTitle("JSON Viewer")

        self.mainLayout = QVBoxLayout(self)  # Set the main layout directly on the widget    
        self.setupDirectoryInput(self.mainLayout)

        self.buttonArrayLayout = QHBoxLayout()
        self.mainLayout.addLayout(self.buttonArrayLayout)

        self.setupTreeView()
        self.setupSwitchViewButton()
        self.setupSaveSettingsButton()
        self.loadInputFiles()

        self.tabWidget = QTabWidget()
        # self.tabWidget.setTabPosition(QTabWidget.West)
        self.setupHistoryList()
        self.setupFileSystemView(self.mainLayout)  # Call the method to set up the file system view
        self.mainLayout.addWidget(self.tabWidget)

        self.setGeometry(100, 100, 1500, 1500)

    def saveCurrentSettings(self):
        # Use the directory chosen by the user as the save directory
        saveDir = self.inputDirEdit.text()
        if not os.path.isdir(saveDir):
            os.makedirs(saveDir)  # Create the directory if it does not exist

        # Save settings.json
        settings_path = os.path.join(saveDir, 'settings.json')
        settings_data = self.settingJSON.toDict()  # Serialize model to a dictionary
        with open(settings_path, 'w') as file:
            json.dump(settings_data, file, indent=4)
        
        # Save input files
        input_files_data = self.inputFilesJSON.toDict()  # Serialize model to a dictionary
        for key, objects in input_files_data.items():
            for obj in objects:
                filename = obj.get("filename", f"default_{key}.json")  # Use a default name if filename is missing
                file_path = os.path.join(saveDir, filename)
                with open(file_path, 'w') as file:
                    json.dump(obj, file, indent=4)  # Save each object literal to its respective file
        
        QMessageBox.information(self, "Save Successful", "Settings and input files have been saved.")

    def setupSaveSettingsButton(self):
        # Setup the save settings button with maximum width and add to buttonArrayLayout
        self.saveSettingsButton = QPushButton("Save Settings")
        self.saveSettingsButton.clicked.connect(self.saveCurrentSettings)
        self.saveSettingsButton.setMaximumWidth(100)
        self.buttonArrayLayout.addWidget(self.saveSettingsButton)       

    def setupDirectoryInput(self, layout):
        """Sets up the directory input field."""
        inputLayout = QHBoxLayout()
        label = QLabel("Input Directory:")
        self.inputDirEdit = ClickableLineEdit()
        self.inputDirEdit.textChanged.connect(self.onDirectoryChanged)

        inputLayout.addWidget(label)
        inputLayout.addWidget(self.inputDirEdit)
        layout.addLayout(inputLayout)

    def showDirectorySummary(self, directoryPath, sortBy='extension'):
        # List to hold file details: name or pattern, extension, and count
        file_details = []

        # List all items in the directory
        for name in os.listdir(directoryPath):
            full_path = os.path.join(directoryPath, name)
            if os.path.isfile(full_path):
                # Use regex to identify the name, numeric part, and extension
                match = re.match(r"(.+?)(\d*)(\.\w+)$", name)
                if match and match.group(2):  # If there is a numeric part
                    # Construct the pattern by excluding the numeric part, include extension for sorting
                    pattern = f"{match.group(1)}*{match.group(3)}"
                    extension = match.group(3)
                    # Find if the pattern already exists in file_details
                    found = False
                    for detail in file_details:
                        if detail[0] == pattern:
                            detail[2] += 1  # Increment count
                            found = True
                            break
                    if not found:
                        file_details.append([pattern, extension, 1])
                else:
                    # For files without a numeric part, use the exact file name
                    extension = match.group(3) if match else ''
                    file_details.append([name, extension, 1])

        # Count directories
        num_dirs = len([name for name in os.listdir(directoryPath) if os.path.isdir(os.path.join(directoryPath, name))])

        # Sorting by the specified criteria
        if sortBy == 'name':
            file_details.sort(key=lambda x: x[0])  # Sort by file name/pattern
        elif sortBy == 'extension':
            file_details.sort(key=lambda x: x[1])  # Sort by extension

        # Prepare the header and table rows
        summary_lines = ["Directory Summary", f"Directory: {directoryPath}", "Pattern/File Name       Count"]
        for detail in file_details:
            summary_lines.append(f"{detail[0].ljust(25)} {detail[2]}")

        summary_lines.append(f"{'Directories:'.ljust(25)} {num_dirs}")
        summary = "\n".join(summary_lines)
        self.fileViewerSummary.setText(summary)


    def setupFileSystemView(self, mainLayout):
        """Sets up the file system view and adds it to the main layout."""
        self.fileSystemModel = QFileSystemModel()
        self.fileSystemModel.setRootPath('')  # Set root path if necessary
        self.fileSystemModel.setFilter(QDir.NoDotAndDotDot | QDir.AllEntries)

        # Create a new QWidget as the container for the layout
        containerWidget = QWidget()
        self.fileViewerLayout = QHBoxLayout(containerWidget)  # Pass the container widget to the layout

        self.fileViewerSummary = QLabel()

        self.fileViewer = QTreeView()
        self.fileViewer.setModel(self.fileSystemModel)
        self.fileViewer.setSortingEnabled(True)
        self.fileViewer.sortByColumn(0, Qt.AscendingOrder)  # Sort by name
        for i in range(1, 4):
            self.fileViewer.hideColumn(i)  # Hide size, type, date modified columns

        # Add the tree views to the layout
        # Assuming self.fileViewerLayout is your QVBoxLayout
        self.fileViewerLayout.addWidget(self.fileViewerSummary)  # Add the summary widget to the layout
        self.fileViewerLayout.addWidget(self.fileViewer)

        # Now add the container widget to the tab
        self.tabWidget.addTab(containerWidget, "File System")
        
    def setupHistoryList(self):
        """Sets up a history list for previously selected directories."""
        self.history_file = "history.json"
        self.historyList = QListWidget()
        self.historyList.itemClicked.connect(self.apply_history)
        # self.layout().addWidget(self.historyList)
        self.load_history()
        self.historyList.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
        self.historyList.setFixedHeight(100)  # Set the fixed height you prefer

        self.saveButton = QPushButton("Save History")
        self.saveButton.clicked.connect(self.save_history)
        self.buttonArrayLayout.addWidget(self.saveButton)
        # left align the button
        self.saveButton.setMaximumWidth(100)
        self.buttonArrayLayout.addStretch()

        self.tabWidget.addTab(self.historyList, "History")
        
    def setupTreeView(self):
        """Sets up the tree views for settings and input files JSON data."""
        self.settingTree = QTreeView()
        self.inputFilesTree = QTreeView()
        
        self.settingJSON = JsonModel()
        self.inputFilesJSON = JsonModel()
        
        self.settingTree.setModel(self.settingJSON)
        self.inputFilesTree.setModel(self.inputFilesJSON)
        
        for tree in (self.settingTree, self.inputFilesTree):
            tree.setAlternatingRowColors(True)
            tree.setTextElideMode(Qt.ElideMiddle)
            tree.setColumnWidth(0, 150)

        self.settingTree.clicked.connect(self.onSettingTreeClicked)

        self.splitter = QSplitter()
        
        self.splitter.addWidget(self.settingTree)
        self.splitter.addWidget(self.inputFilesTree)
        self.layout().addWidget(self.splitter)

        total_width = self.splitter.width()        
        setting_tree_width = int(total_width * 4 / (4 + 6))  # 3 parts of the total width
        input_files_tree_width = int(total_width * 6 / (4 + 6))  # 7 parts of the total width
        self.splitter.setSizes([setting_tree_width, input_files_tree_width])
        self.splitter.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)  # Ensure splitter expands
        self.mainLayout.addWidget(self.splitter)  # Add splitter to the main layout

    def setupSwitchViewButton(self):
        # Setup the switch view button with maximum width and add to buttonArrayLayout
        self.switchViewButton = QPushButton("Switch")
        self.switchViewButton.clicked.connect(self.switchTreeViewOrientation)
        self.switchViewButton.setMaximumWidth(100)
        self.buttonArrayLayout.addWidget(self.switchViewButton)

    def onSettingTreeClicked(self, index):
        itemData = index.data(Qt.DisplayRole)  # The visible text of the item
        # Attempt to find the item in the inputFilesJSON model
        foundIndex = self.inputFilesJSON.findInInputFilesJSON(itemData)
        if foundIndex.isValid():
            self.inputFilesTree.setCurrentIndex(foundIndex)
            self.inputFilesTree.scrollTo(foundIndex, QTreeView.PositionAtTop)
        else:
            if os.path.isdir(itemData):
                self.fileViewer.setRootIndex(self.fileSystemModel.setRootPath(itemData))
                self.showDirectorySummary(itemData)  # Display summary of the directory
            else:
                pass  # Do nothing or show a message

    def switchTreeViewOrientation(self):
        currentOrientation = self.splitter.orientation()
        if currentOrientation == Qt.Horizontal:
            self.splitter.setOrientation(Qt.Vertical)
        else:
            self.splitter.setOrientation(Qt.Horizontal)

    def save_history(self):
        input_dir = self.inputDirEdit.text()
        if not os.path.isdir(input_dir):
            QMessageBox.warning(self, "Error", "Invalid directory.")
            return
        try:
            history = []
            if os.path.exists(self.history_file):
                with open(self.history_file, "r") as file:
                    history = json.load(file)
            # Check if the directory is already in history
            if not any(item['path'] == input_dir for item in history):
                history.append({
                    "path": input_dir,
                    "date_added": str(datetime.date.today())
                })
                with open(self.history_file, "w") as file:
                    json.dump(history, file, indent=4)
                self.update_history_list()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save history: {str(e)}")

    def update_history_list(self):
        self.historyList.clear()
        if os.path.exists(self.history_file):
            with open(self.history_file, "r") as file:
                history = json.load(file)
                for entry in history:
                    self.historyList.addItem(f"{entry['path']} (Added on: {entry['date_added']})")

    def apply_history(self, item):
        # Parse the directory path from the list item's text
        path = item.text().split(" (Added on:")[0]
        self.inputDirEdit.setText(path)

    def load_history(self):
        self.update_history_list()

    def onDirectoryChanged(self, newDirPath):
        if newDirPath:  # Check if the new directory path is not empty
            self.loadInputFiles(newDirPath)

    def loadInputFiles(self, dir_path = './input_files/Tanizwa1996_H0d05_L1d8_piston/'):
        settings_path = f"{dir_path}/setting.json"
        print(settings_path)
        settings = {}
        input_files = []
        try:
            with open(settings_path, 'r') as file:
                settings = json.load(file)
                input_files = settings.get("input_files", [])
        except FileNotFoundError:
            QMessageBox.critical(self, "Error", "Settings file not found.")
            return

        print(input_files)

        type_documents_filename = {}
        for file_name in input_files:
            try:
                with open(f"{dir_path}/{file_name}", 'r') as file:
                    document = json.load(file)
                    doc_type = document.get('type')
                    if doc_type:
                        type_documents_filename.setdefault(doc_type, []).append({"filename": file_name, **document})
            except FileNotFoundError:
                QMessageBox.critical(self, "Error", f"File {file_name} not found.")
                continue

        self.inputFilesJSON.load(type_documents_filename)
        self.inputFilesTree.expandAll()

        self.settingJSON.load(settings)
        self.settingTree.expandAll()

    # Inside your class
    def selectInputDirectory(self):
        dir = QFileDialog.getExistingDirectory(self, "Select Input Directory")
        if dir:  # Check if a directory was selected
            self.inputFileEdit.setText(dir)

    def selectOutputDirectory(self):
        dir = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if dir:  # Check if a directory was selected
            self.outputDirEdit.setText(dir)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()

    window.show()
    sys.exit(app.exec())
