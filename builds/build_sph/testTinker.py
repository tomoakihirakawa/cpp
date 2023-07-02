import subprocess
import tkinter as tk
from tkinter import ttk


def run_command():
    # Run your command
    process = subprocess.Popen(
        ["./main"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    # Update the text widget with the command output
    text_widget.insert(tk.END, stdout)


root = tk.Tk()

# Create a button to run the command
button = ttk.Button(root, text="Run Command", command=run_command)
button.pack()

# Create a text widget to show the command output
text_widget = tk.Text(root, width=50, height=10)
text_widget.pack()

root.mainloop()
