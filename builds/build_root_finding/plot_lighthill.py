import os
import matplotlib.pyplot as plt
import imageio.v2 as imageio
import numpy as np

plt.rcParams["font.family"] = "Times New Roman"


def plot_data(file_name, analitical_file_name):
    i_data, x_data, y_data = [], [], []

    # Read data from file
    with open(file_name, 'r') as f:
        for line in f:
            index, x, y, q = map(float, line.strip().split())
            i_data.append(index)
            x_data.append(x)
            y_data.append(y)

    # Read analytical data
    x_analitical, y_analitical = [], []
    with open(analitical_file_name, 'r') as f:
        for line in f:
            x, y = map(float, line.strip().split())
            x_analitical.append(x)
            y_analitical.append(y)

    # Create plot
    fig, ax = plt.subplots(figsize=(8, 3))

    # Set range
    ax.set_xlim(0., 0.8)
    ax.set_ylim(-0.1, 0.1)

    # Plot point with larger size
    ax.scatter(x_data[1:-1], y_data[1:-1], s=50, label='node')

    for i, txt in enumerate(y_data[1:-1]):
        ax.text(x_data[i+1], y_data[i+1],
                str(int(i_data[i+1])), fontsize=12, color='red')

    # Aspect ratio
    ax.set_aspect('equal', adjustable='box')

    # Grid
    ax.grid()

    # Plot line
    ax.plot(x_data, y_data, label='Numerical Solution')

    # Plot analytical solution
    ax.plot(x_analitical, y_analitical, label='Analytical Solution')

    # Labels and title with larger font size and times
    ax.set_xlabel('x', fontsize=12, fontname='Times New Roman')
    ax.set_ylabel('y', fontsize=12, fontname='Times New Roman')
    ax.set_title(file_name, fontsize=14, fontname='Times New Roman')

    # Legend
    ax.legend(loc='upper right')

    # Save the figure
    plt.savefig(file_name + ".png")

    # Close the figure
    plt.close()


# Directory containing data files
dir_name = "./output_lighthill/"

# List of data files
files = [f for f in os.listdir(dir_name) if f.endswith(
    '.txt') and not f.endswith('_analitical.txt')]

# Sort files to ensure they're processed in order
files.sort(key=lambda x: int(x.replace('lighthill', '').replace('.txt', '')))

# List to hold images for GIF
images = []

for file_name in files:
    analitical_file_name = file_name.replace('.txt', '_analitical.txt')
    plot_data(dir_name + file_name, dir_name + analitical_file_name)
    images.append(imageio.imread(dir_name + file_name + ".png"))

imageio.mimsave('sample.gif', images, loop=100)
