import os
import matplotlib.pyplot as plt
import imageio.v2 as imageio  # This is the corrected line

# Function to read data from a file and plot it

plt.rcParams["font.family"] = "Times New Roman"

# Function to read data from a file and plot it


def plot_data(file_name):
    x_data, y_data = [], []

    # Read data from file
    with open(file_name, 'r') as f:
        for line in f:
            x, y, q = map(float, line.strip().split())
            x_data.append(x)
            y_data.append(y)

    # Create plot
    fig, ax = plt.subplots(figsize=(8, 3))

    # Set range
    ax.set_xlim(0., 0.8)
    ax.set_ylim(-0.1, 0.1)

    # Plot point with larger size
    ax.scatter(x_data, y_data, s=50, label='node')

    # Aspect ratio
    ax.set_aspect('equal', adjustable='box')

    # Grid
    ax.grid()

    # Plot line
    ax.plot(x_data, y_data)

    # Labels and title with larger font size and times
    ax.set_xlabel('x', fontsize=12, fontname='Times New Roman')
    ax.set_ylabel('y', fontsize=12, fontname='Times New Roman')
    ax.set_title(file_name, fontsize=14, fontname='Times New Roman')

    # Legend
    ax.legend()

    # Save the figure
    plt.savefig(file_name + ".png")

    # Close the figure
    plt.close()


# Directory containing data files
dir_name = "./output_lighthill/"

# List of data files
files = [f for f in os.listdir(dir_name) if f.endswith('.txt')]

# Sort files to ensure they're processed in order
files.sort(key=lambda x: int(x.replace('lighthill', '').replace('.txt', '')))

# List to hold images for GIF
images = []

for file_name in files:
    plot_data(dir_name + file_name)
    # This is the corrected line
    images.append(imageio.imread(dir_name + file_name + ".png"))

# Save as a GIF
imageio.mimsave(dir_name + 'robot_movement.gif', images)
