'''

```shell
pip install -U PySide6
```
する

'''

import sys
import os
import json
from PySide6.QtWidgets import (
    QApplication, QWidget, QLineEdit, QFormLayout, QPushButton, QVBoxLayout, QTextBrowser, QListWidget, QMessageBox
)
from PySide6.QtCore import QProcess
from PySide6.QtWidgets import QMenu
from PySide6.QtCore import QPoint
from PySide6.QtCore import Qt

def ansi_to_html(text):
    replacements = {

        "\033[0;31m": "<span style='color:red;'>",  # Red
        # Font colors
        "\033[30m": "<span style='color:black;'>",  # Black
        "\033[31m": "<span style='color:red;'>",  # Red
        "\033[32m": "<span style='color:green;'>",  # Green
        "\033[33m": "<span style='color:yellow;'>",  # Yellow
        "\033[34m": "<span style='color:blue;'>",  # Blue
        "\033[35m": "<span style='color:magenta;'>",  # Magenta
        "\033[36m": "<span style='color:cyan;'>",  # Cyan
        "\033[37m": "<span style='color:white;'>",  # White
        
        # Bright colors
        "\033[1;30m": "<span style='color:black; font-weight:bold;'>",  # Bright Black (Bold)
        "\033[1;31m": "<span style='color:red; font-weight:bold;'>",  # Bright Red (Bold)
        "\033[1;32m": "<span style='color:green; font-weight:bold;'>",  # Bright Green (Bold)
        "\033[1;33m": "<span style='color:yellow; font-weight:bold;'>",  # Bright Yellow (Bold)
        "\033[1;34m": "<span style='color:blue; font-weight:bold;'>",  # Bright Blue (Bold)
        "\033[1;35m": "<span style='color:magenta; font-weight:bold;'>",  # Bright Magenta (Bold)
        "\033[1;36m": "<span style='color:cyan; font-weight:bold;'>",  # Bright Cyan (Bold)
        "\033[1;37m": "<span style='color:white; font-weight:bold;'>",  # Bright White (Bold)
        
        # Styles
        "\033[0m": "</span>",  # Reset
        "\033[1m": "<span style='font-weight:bold;'>",  # Bold
        "\033[3m": "<span style='font-style:italic;'>",  # Italic
        "\033[4m": "<span style='text-decoration:underline;'>",  # Underline
        
        # Background Colors
        "\033[40m": "<span style='background-color:black;'>",  # Background Black
        "\033[41m": "<span style='background-color:red;'>",  # Background Red
        "\033[42m": "<span style='background-color:green;'>",  # Background Green
        "\033[43m": "<span style='background-color:yellow;'>",  # Background Yellow
        "\033[44m": "<span style='background-color:blue;'>",  # Background Blue
        "\033[45m": "<span style='background-color:magenta;'>",  # Background Magenta
        "\033[46m": "<span style='background-color:cyan;'>",  # Background Cyan
        "\033[47m": "<span style='background-color:white;'>",  # Background White
    }
    # Replace ANSI codes with HTML equivalent
    for ansi, html in replacements.items():
        text = text.replace(ansi, html)
    # Close all open HTML tags at the end
    text += "</span>" * text.count("<span")
    return text

class RemeshWidget(QWidget):
    def __init__(self):
        super().__init__()

        self.setLayout(QFormLayout())
        self.inputFileEdit = QLineEdit()
        self.outputDirEdit = QLineEdit()
        self.fileNameEdit = QLineEdit()
        self.remeshNumberEdit = QLineEdit()
        self.remeshButton = QPushButton("Remesh")
        self.saveButton = QPushButton("Save History")  # New save button
        self.consoleOutput = QTextBrowser()
        self.historyList = QListWidget()  # Widget to display history

        self.inputFileEdit.setFixedSize(700, 30)
        self.outputDirEdit.setFixedSize(700, 30)
        self.consoleOutput.setFixedSize(700, 300)
        self.remeshButton.setFixedSize(100, 30)
        self.saveButton.setFixedSize(100, 30)
        self.historyList.setFixedSize(700, 200)  # Adjust size as needed
        self.consoleOutput.setOpenExternalLinks(True)  # Allow HTML
        self.consoleOutput.setAcceptRichText(True)  # Interpret as rich text

        self.layout().addRow("Input File:", self.inputFileEdit)
        self.layout().addRow("Output Directory:", self.outputDirEdit)
        self.layout().addRow("File Name:", self.fileNameEdit)
        self.layout().addRow("Remesh Number:", self.remeshNumberEdit)
        self.layout().addRow(self.remeshButton)
        self.layout().addRow(self.saveButton)  # Add the save button to the layout
        self.layout().addRow("Console Output:", self.consoleOutput)
        self.layout().addRow("History:", self.historyList)  # Add the history list to the layout

        self.stopButton = QPushButton("Stop Process")
        self.layout().addRow(self.stopButton)
        self.stopButton.clicked.connect(self.stop_process)

        self.remeshButton.clicked.connect(self.run_remesh)
        self.saveButton.clicked.connect(self.save_history)  # Connect save button to save_history function
        self.process = QProcess(self)
        self.process.readyReadStandardOutput.connect(self.handle_stdout)
        self.process.readyReadStandardError.connect(self.handle_stderr)
        self.historyList.itemClicked.connect(self.apply_history)
        self.historyList.setContextMenuPolicy(Qt.CustomContextMenu)
        self.historyList.customContextMenuRequested.connect(self.on_context_menu)

        self.history_file = "remesh_history.json"
        self.load_history()  # Load history at startup

    def stop_process(self):
        if self.process and self.process.state() != QProcess.NotRunning:
            self.process.terminate()  # Gracefully terminate the process
            # Optionally, wait for a specific time and then forcefully kill if still running
            if not self.process.waitForFinished(3000):  # Wait for 3 seconds
                self.process.kill()  # Forcefully kill the process if it didn't terminate
                self.consoleOutput.append("Process was forcefully stopped.")
            else:
                self.consoleOutput.append("Process stopped successfully.")
        else:
            self.consoleOutput.append("No process is currently running.")

    def run_remesh(self):
        input_file = self.inputFileEdit.text()
        output_dir = self.outputDirEdit.text()
        file_name = self.fileNameEdit.text()
        remesh_number = self.remeshNumberEdit.text()
        if not all([input_file, output_dir, file_name, remesh_number]):
            self.consoleOutput.append("Error: All fields must be filled out.")
            return
        command = f"./remesh {input_file} {output_dir} {file_name} {remesh_number}"
        self.process.start("zsh", ["-c", command])
        # Removed automatic saving to history

    def save_history(self):
        # Modified to save current GUI input fields instead of parameters
        input_file = self.inputFileEdit.text()
        output_dir = self.outputDirEdit.text()
        file_name = self.fileNameEdit.text()
        remesh_number = self.remeshNumberEdit.text()
        try:
            history = []
            if os.path.exists(self.history_file):
                with open(self.history_file, "r") as file:
                    history = json.load(file)

            history.append({
                "input_file": input_file,
                "output_dir": output_dir,
                "file_name": file_name,
                "remesh_number": remesh_number
            })

            with open(self.history_file, "w") as file:
                json.dump(history, file, indent=4)

            self.update_history_list()  # Update the history list widget with new entry
        except Exception as e:
            self.consoleOutput.append(f"Failed to save history: {str(e)}")

    def update_history_list(self):
        # Clear the list and repopulate with updated history
        self.historyList.clear()
        if os.path.exists(self.history_file):
            with open(self.history_file, "r") as file:
                history = json.load(file)
            for operation in history:
                self.historyList.addItem(f"{operation['input_file']} -> {operation['output_dir']}")

    def load_history(self):
        self.update_history_list()  # Use the new method to load history into the list widget

    # Example logging for debugging
    def handle_stdout(self):
        data = self.process.readAllStandardOutput().data().decode()
        html_text = ansi_to_html(data)
        # Use insertHtml to add formatted HTML directly without replacing existing content
        self.consoleOutput.insertHtml(html_text)
        # Ensure there's a new line after the inserted HTML content
        self.consoleOutput.append('')

    def handle_stderr(self):
        error = self.process.readAllStandardError().data().decode()
        self.consoleOutput.append(error)  # Append as HTML

    def apply_history(self, item):
        with open(self.history_file, 'r') as file:
            history = json.load(file)
        
        for entry in history:
            formatted_entry = f"{entry['input_file']} -> {entry['output_dir']}"
            if formatted_entry == item.text():
                self.inputFileEdit.setText(entry['input_file'])
                self.outputDirEdit.setText(entry['output_dir'])
                self.fileNameEdit.setText(entry['file_name'])
                self.remeshNumberEdit.setText(entry['remesh_number'])
                break

    def on_context_menu(self, position):
        contextMenu = QMenu(self)
        deleteAction = contextMenu.addAction("Delete")
        action = contextMenu.exec(self.historyList.mapToGlobal(position))
        if action == deleteAction:
            self.delete_selected_history()

    def delete_selected_history(self):
        currentItem = self.historyList.currentItem()
        if currentItem is None:
            self.consoleOutput.append("No history selected for deletion.")
            return        
        # Confirm before deletion
        reply = QMessageBox.question(self, 'Delete History', 'Are you sure you want to delete this history?',
                                    QMessageBox.Yes | QMessageBox.No, QMessageBox.No)

        if reply == QMessageBox.Yes:
            # Find and remove the selected entry from the JSON file
            with open(self.history_file, 'r') as file:
                history = json.load(file)
            
            # Assuming each list item corresponds exactly to a history entry
            entry_to_delete = currentItem.text()
            for i, entry in enumerate(history):
                formatted_entry = f"{entry['input_file']} -> {entry['output_dir']}"
                if formatted_entry == entry_to_delete:
                    del history[i]
                    break
            
            # Save the updated history back to the JSON file
            with open(self.history_file, 'w') as file:
                json.dump(history, file, indent=4)
            
            # Remove the entry from the list widget
            self.historyList.takeItem(self.historyList.row(currentItem))
            self.consoleOutput.append("Selected history deleted successfully.")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    widget = RemeshWidget()
    widget.show()
    sys.exit(app.exec())
