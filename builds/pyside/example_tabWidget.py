import sys
from PySide6.QtWidgets import QApplication, QMainWindow, QTabWidget, QLabel, QVBoxLayout, QWidget

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Tab Widget Example")
        
        # Create the tab widget
        self.tabWidget = QTabWidget()
        
        # Add tabs
        self.addTab("Tab 1", "Content of Tab 1")
        self.addTab("Tab 2", "Content of Tab 2")
        self.addTab("Tab 3", "Content of Tab 3")
        
        # Set the central widget
        self.setCentralWidget(self.tabWidget)
        
    def addTab(self, title, content):
        # Create a widget and layout for the tab content
        tabContentWidget = QWidget()
        layout = QVBoxLayout()
        tabContentWidget.setLayout(layout)
        
        # Add a label as the content of the tab
        label = QLabel(content)
        layout.addWidget(label)
        
        # Add the tab to the tab widget
        self.tabWidget.addTab(tabContentWidget, title)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    
    mainWindow = MainWindow()
    mainWindow.show()
    
    sys.exit(app.exec())
