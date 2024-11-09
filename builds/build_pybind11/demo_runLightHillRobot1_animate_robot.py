
'''DOC_EXTRACT 0_2_0_result

## python内で共有ライブラリを使う

### アニメーションgifファイルを作成しロボットの動きを可視化する

\ref{PYBIND11:HOW_TO_IMPORT}{このように}`import`して利用できる．
cppと同じように\ref{PYBIND11:GENERATE_OBJECT}{`robot`オブジェクトを作成}．

出力結果

|例：水族館の魚|出力結果|
|:---:|:---:|
| <img src="sample_aquarium.gif"  width="80%" height="80%"> | ![sample.gif](sample.gif) |

'''

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
# \label{PYBIND11:HOW_TO_IMPORT}
import LighthillRobot as LHR

# ---------------------------------- ロボットの設定 --------------------------------- #

L = 0.25
period = 1.
w = 2.*math.pi/period
k = 2.*math.pi / L
c1 = 0.02
c2 = 0.01
n = 2

# \label{PYBIND11:GENERATE_OBJECT}
robot = LHR.LighthillRobot(L, w, k, c1, c2, n)

# -------------------------------- gifファイルの作成 -------------------------------- #
fig, ax = plt.subplots()

# set aspect ratio as real size
ax.set_aspect('equal')

line, = ax.plot([], [], lw=2, marker='.')

def init():
    ax.set_xlim(0, L)
    ax.set_ylim(-L/3., L/3.)
    return line,

text_handle = None  # Initialize text handle outside update function

def update(t):
    global text_handle  # Declare as global to modify

    # robot.c1 = 0.2*math.sin(t)
    # robot.c2 = 0.2*math.sin(t)
    angles = robot.getAngles(t)
    c = (math.tanh((t-1.)*math.pi)+1.)/2.
    angles[0] *= c
    angles[1] *= c
    angles[2] *= c
    positions = robot.anglesToX(angles)
    positions_x = [pos[0] for pos in positions]
    positions_y = [pos[1] for pos in positions]

    line.set_data(positions_x, positions_y)
    # You can adjust the position
    if text_handle is not None:
        text_handle.remove()
    text_handle = ax.text(
        0.1, 0.9, f'Time: {t:.2f}', transform=ax.transAxes)

    return line,

# Define time array
T_max = 10.0  # Maximum time
dt = 0.05  # Time step size
t_values = np.arange(0, T_max, dt)  # Array of time values

ani = animation.FuncAnimation(
    fig, update, frames=t_values, init_func=init, blit=True, interval=1000/50)
ani.save("sample.gif", writer='pillow')

plt.show()

