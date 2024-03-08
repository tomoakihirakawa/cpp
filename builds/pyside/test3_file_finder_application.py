import fnmatch
import sys
import os
from PySide6.QtWidgets import QApplication, QMainWindow, QWidget, QVBoxLayout, QPushButton, QLineEdit, QListWidget, QFileDialog

class FileFinder(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("File Finder")
        
        # Main layout
        layout = QVBoxLayout()
        
        # Widgets
        self.pathEdit = QLineEdit()
        self.pathEdit.setPlaceholderText("Choose directory...")
        self.browseButton = QPushButton("Browse")
        self.patternEdit = QLineEdit("*.txt")
        self.findButton = QPushButton("Find")
        self.resultList = QListWidget()
        
        # Adding widgets to layout
        layout.addWidget(self.pathEdit)
        layout.addWidget(self.browseButton)
        layout.addWidget(self.patternEdit)
        layout.addWidget(self.findButton)
        layout.addWidget(self.resultList)
        
        # Setting the main widget
        centralWidget = QWidget()
        centralWidget.setLayout(layout)
        self.setCentralWidget(centralWidget)
        
        # Connect signals
        self.browseButton.clicked.connect(self.browseDirectory)
        self.findButton.clicked.connect(self.findFiles)

    def browseDirectory(self):
        directory = QFileDialog.getExistingDirectory(self, "Select Directory")
        if directory:
            self.pathEdit.setText(directory)
    
    def findFiles(self):
        self.resultList.clear()
        directory = self.pathEdit.text()
        pattern = self.patternEdit.text()
        if directory and pattern:
            for root, dirs, files in os.walk(directory):
                for file in files:
                    if self.patternMatch(file, pattern):
                        self.resultList.addItem(os.path.join(root, file))
    
    def patternMatch(self, filename, pattern):
        # Simple pattern matching for demonstration, replace with more complex logic if needed
        return fnmatch.fnmatch(filename, pattern)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = FileFinder()
    window.show()
    sys.exit(app.exec())
