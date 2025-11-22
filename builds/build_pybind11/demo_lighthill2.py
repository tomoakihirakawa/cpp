import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

L = 1.0
w = 2.0 * np.pi
k = 2.0 * np.pi
c1 = 0.2
c2 = 0.2


def yLH(x, t):
    return (c1 * x / L + c2 * (x / L)**2) * np.sin(k * (x / L) - w * t)


# Create a new figure and add an axis
fig, ax = plt.subplots()

# Create an array of x values from 0 to L
x = np.linspace(0, L, 1000)

# Create a line object. We will update this line in the animate function
line, = ax.plot(x, yLH(x, 0))

# Set up the axis limits
ax.set_xlim(-L*0.5, L+L*0.5)
ax.set_ylim(-1, 1)

# Define the animation function


def animate(i):
    line.set_ydata(yLH(x, i * 0.01))  # update the data
    return line,


# Create the animation
ani = FuncAnimation(fig, animate, frames=1000, interval=10, blit=True)

# Display the animation
plt.show()
