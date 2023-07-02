# pip install pyvista
#
import os
from tkinter import ttk, filedialog
import tkinter as tk
import json
import pyvista as pv


root = tk.Tk()  # Initialize tkinter root widget
root.title('My Application')
root.geometry('800x600')  # Set appropriate size for your window

# String variables for entry fields
name_var = tk.StringVar(root)
type_var = tk.StringVar(root)
objfile_var = tk.StringVar(root)

# Define object types and properties
object_types = ["Fluid", "RigidBody", "SoftBody"]
properties = {
    "Fluid": ["name", "type", "objfile", "output_vtu_file_name", "output_pvd_file_name"],
    "RigidBody": ["name", "type", "objfile", "output_vtu_file_name", "output_pvd_file_name"],
    "SoftBody": ["name", "type", "objfile", "output_vtu_file_name", "output_pvd_file_name"]
}

# Define column titles for tree view
titles = ("Name", "Type", "Obj file")

# Initialize Treeview widget
objects_tree = ttk.Treeview(root, columns=titles, show='headings')
for title in titles:
    objects_tree.heading(title, text=title)
objects_tree.column("Name", width=100)
objects_tree.column("Type", width=100)
objects_tree.column("Obj file", width=600)
objects_tree.pack(padx=30, pady=30)  # Position the Treeview widget

objects_dict = {}  # Dictionary to hold objects

# Frame for 'add object' interface
add_frame = ttk.Frame(root, padding="30 30 30 30")
name_ttk_label = ttk.Label(add_frame, text=titles[0])
type_ttk_label = ttk.Label(add_frame, text=titles[1])
file_ttk_label = ttk.Label(add_frame, text=titles[2])
name_entry = ttk.Entry(add_frame, textvariable=name_var)
type_entry = ttk.Combobox(add_frame, textvariable=type_var,
                          values=object_types, state='readonly')
file_entry = ttk.Entry(add_frame, textvariable=objfile_var)
file_entry.bind("<Button-1>", lambda event: objfile_var.set(filedialog.askopenfilename(
    initialdir=os.getcwd(), filetypes=(("Obj files", "*.obj"), ("All files", "*.*")))))

# Positioning 'add object' interface
name_ttk_label.grid(row=0, column=0)
type_ttk_label.grid(row=0, column=1)
file_ttk_label.grid(row=0, column=2)
name_entry.grid(row=1, column=0)
type_entry.grid(row=1, column=1)
file_entry.grid(row=1, column=2)
add_frame.pack()

# Function to add new object


def add_object():
    # Fetch details from the user
    name = name_var.get()
    obj_type = type_var.get()
    obj_file = objfile_var.get()

    # Check if all required fields have values
    if not all([name, obj_type, obj_file]):
        print("Please fill all fields.")
        return

    # Create new object and add to dictionary and treeview
    new_obj = {prop: name if prop == "name" else (
        obj_file if prop == "objfile" else "") for prop in properties[obj_type]}
    objects_dict[name] = new_obj
    objects_tree.insert("", tk.END, values=(name, obj_type, obj_file))


# Button to add new object
add_button = ttk.Button(add_frame, text="add", command=add_object)
add_button.grid(row=1, column=3)

# Frame for object properties
properties_frame = tk.Frame(root)
properties_frame.pack_forget()  # Initially hidden

# Function to add a property field


def add_property_entry(i, prop, obj):
    tk.Label(properties_frame, text=prop).grid(row=i, column=0)
    var = tk.StringVar(root)
    var.set(obj.get(prop, ''))
    entry = tk.Entry(properties_frame, textvariable=var, width=100)
    entry.grid(row=i, column=1)
    # Update object when property changes
    var.trace("w", lambda *args: obj.update({prop: var.get()}))

# Function to load properties of an object when selected in the treeview


def view_3d_model(path_to_file):
    # Load the STL file using pyvista
    mesh = pv.read(path_to_file)

    # Create a Plotter object
    plotter = pv.Plotter(notebook=False)

    # Add the mesh to the plotter
    plotter.add_mesh(mesh)

    # Render the plotter
    plotter.show()


def load_properties(event):
    selected_item = objects_tree.selection()[0]
    selected_name = objects_tree.item(selected_item)['values'][0]
    obj = objects_dict[selected_name]

    # Clear current properties
    for widget in properties_frame.winfo_children():
        widget.destroy()

    # Add property fields for each property
    for i, prop in enumerate(obj):
        add_property_entry(i, prop, obj)

    obj_file = objects_dict[selected_name]['objfile']
    view_3d_model(obj_file)

    properties_frame.pack(padx=30, pady=30)
    root.update_idletasks()  # Force redraw

# Function to load settings from a json file


def load_settings():
    settings_file = filedialog.askopenfilename(initialdir=os.getcwd(
    ), filetypes=(("JSON files", "*.json"), ("All files", "*.*")))

    if settings_file:
        settings_dir = os.path.dirname(settings_file)

        with open(settings_file, 'r') as f:
            settings = json.load(f)

        # Load each object from a separate file
        for json_file in settings.get('input_files', []):
            json_file_path = os.path.join(settings_dir, json_file)

            with open(json_file_path, 'r') as f:
                obj = json.load(f)

            # Add object to dictionary and treeview
            if 'name' in obj and 'type' in obj and 'objfile' in obj:
                name = obj['name']
                objects_dict[name] = obj
                objects_tree.insert("", tk.END, values=(
                    obj['name'], obj['type'], obj['objfile']))
            else:
                print(
                    f'Object in {json_file} missing "name", "type", or "objfile"')


# Button to load settings
load_button = ttk.Button(add_frame, text="load", command=load_settings)
load_button.grid(row=1, column=4)

# Bind selection event of treeview to load properties function
objects_tree.bind("<<TreeviewSelect>>", load_properties)

root.mainloop()  # Start tkinter event loop
