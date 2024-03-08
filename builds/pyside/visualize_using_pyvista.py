import sys
import pyvista as pv
from PySide6.QtWidgets import QApplication, QWidget, QVBoxLayout, QPushButton, QFileDialog, QSlider, QLabel
from PySide6.QtCore import Qt

class MainWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("OBJ File Visualizer")
        self.layout = QVBoxLayout()
        self.button_open = QPushButton("Open OBJ File")
        self.button_add = QPushButton("Add Another OBJ File")
        self.button_visualize = QPushButton("Visualize")
        self.button_visualize.setEnabled(False)  # Disable until file is selected
        self.button_add.setEnabled(False)  # Disable until initial visualization

        self.opacity_label = QLabel("Adjust Opacity of Selected Mesh:")
        self.opacity_slider = QSlider(Qt.Horizontal)
        self.opacity_slider.setMinimum(0)
        self.opacity_slider.setMaximum(100)
        self.opacity_slider.setValue(100)  # Start at 100% opacity
        self.opacity_slider.setEnabled(False)  # Disable until a mesh is picked
        self.opacity_slider.valueChanged.connect(self.adjust_opacity)

        self.button_erase = QPushButton("Erase Selected Mesh")
        self.button_erase.clicked.connect(self.erase_selected_mesh)
        self.button_erase.setEnabled(False)  # Disable until a mesh is picked

        self.button_open.clicked.connect(self.open_file)
        self.button_visualize.clicked.connect(self.visualize)
        self.button_add.clicked.connect(self.add_another_obj)

        self.layout.addWidget(self.button_open)
        self.layout.addWidget(self.button_visualize)
        self.layout.addWidget(self.button_add)
        self.layout.addWidget(self.opacity_label)
        self.layout.addWidget(self.opacity_slider)
        self.layout.addWidget(self.button_erase)
        self.setLayout(self.layout)

        self.obj_file_path = None
        self.plotter = None  # Initialize the plotter as None
        self.selected_mesh = None  # Track the selected mesh

    def open_file(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Open OBJ File", "", "OBJ Files (*.obj)")
        if file_path:
            self.obj_file_path = file_path
            self.button_visualize.setEnabled(True)

    def visualize(self):
        if self.obj_file_path:
            mesh = pv.read(self.obj_file_path)
            if not self.plotter:
                self.plotter = pv.Plotter()
                self.plotter.enable_cell_picking(callback=self.on_pick, show_message=True, style='surface')
                self.plotter.add_mesh(mesh, show_edges=True)
                self.plotter.show(interactive_update=True, auto_close=False)
                self.button_add.setEnabled(True)  # Enable adding more OBJ files after initial visualization
            else:
                self.plotter.add_mesh(mesh, show_edges=True)
                self.plotter.render()

    def add_another_obj(self):
        # Reuse the open file dialog to select another OBJ file
        self.open_file()
        # Automatically visualize the newly selected OBJ without pressing the Visualize button again
        self.visualize()

    # def on_pick(self, picker):
    #     self.selected_mesh = picker.GetDataSet()
    #     self.opacity_slider.setEnabled(True)
    #     self.button_erase.setEnabled(True)

    def on_pick(self, pick_result):
        # Assuming pick_result gives you information about the pick event
        # This doesn't directly solve selecting and modifying individual meshes
        print("Picked at position:", pick_result.position)


    def adjust_opacity(self, value):
        if self.selected_mesh:
            opacity = value / 100.0
            self.plotter.update_scalars(opacity, mesh=self.selected_mesh, render=True)

    def erase_selected_mesh(self):
        if self.selected_mesh and self.plotter:
            self.plotter.remove_actor(self.selected_mesh)
            self.selected_mesh = None
            self.opacity_slider.setEnabled(False)
            self.button_erase.setEnabled(False)
            self.plotter.render()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
