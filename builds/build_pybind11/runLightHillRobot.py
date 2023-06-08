'''DOC_EXTRACT 

出力結果

|水族館の魚|ロボット|
|:---:|:---:|
| ![sample1.gif](sample1.gif) | ![sample.gif](sample.gif) |

'''

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
# \label{PYBIND11:HOW_TO_IMPORT}
import LightHillRobot as LHR

L = .71
w = 2.*math.pi
k = 2.*math.pi
c1 = 0.2
c2 = 0.2
n = 10

# Create a LightHillRobot object
robot = LHR.LightHillRobot(L, w, k, c1, c2, n)

fig, ax = plt.subplots()

# Initialize a line object for the plot
ax.set_aspect('equal', 'box')

line, = ax.plot([], [], lw=2, marker='.')


def init():
    ax.set_xlim(0, L)
    ax.set_ylim(-L/2., L/2.)
    return line,


def update(t):
    # robot.c1 = 0.2*math.sin(t)
    # robot.c2 = 0.2*math.sin(t)
    angles = robot.getAngles(t)
    positions = robot.anglesToX(angles)
    positions_x = [pos[0] for pos in positions]
    positions_y = [pos[1] for pos in positions]

    line.set_data(positions_x, positions_y)
    return line,


# Define time array
T_max = 10.0  # Maximum time
dt = 0.01  # Time step size
t_values = np.arange(0, T_max, dt)  # Array of time values

ani = animation.FuncAnimation(
    fig, update, frames=t_values, init_func=init, blit=True, interval=1000/50)

ani.save("sample.gif", writer='pillow')

plt.show()
