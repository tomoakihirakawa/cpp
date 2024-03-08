import sys
from PySide6.QtWidgets import (
    QApplication,
    QWidget,
    QPushButton,
    QVBoxLayout,
    QFileDialog,
)

# import dir paraview in this directory

import simple as pvs


# ------  Load Initial OBJ --------
# file_path = "bunny.obj"
# reader = pvs.OBJReader(FileName=file_path)
# display = pvs.Show(reader)


# # ------  Button Functionality --------
# def add_obj():
#     file_name, _ = QFileDialog.getOpenFileName(
#         None, "Select OBJ File", "", "OBJ Files (*.obj)"
#     )
#     if file_name:
#         new_reader = pvs.OBJReader(FileName=file_name)
#         pvs.Show(new_reader)
#         pvs.Render() 


# # ------  GUI Setup --------
# app = QQApplication(sys.argv)
# widget = QWidget()
# layout = QVBoxLayout(widget)

# add_button = QPushButton("Add OBJ")
# add_button.clicked.connect(add_obj)
# layout.addWidget(add_button)

# widget.show()

# ----- ParaView Control -----
# pvs.ResetCamera()
# pvs.Render()

# app.exec()
