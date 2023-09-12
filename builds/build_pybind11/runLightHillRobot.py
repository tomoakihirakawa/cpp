
'''DOC_EXTRACT 0_2_0_result

## python内で共有ライブラリを使う

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

L = 1.
period = 1.
w = 2.*math.pi/period
k = 2.*math.pi
c1 = 0.05
c2 = 0.05
n = 2

# \label{PYBIND11:GENERATE_OBJECT}
robot = LHR.LighthillRobot(L, w, k, c1, c2, n)

# -------------------------------- gifファイルの作成 -------------------------------- #
if True:
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

# -------------------------------- datファイルの作成 -------------------------------- #

'''DOC_EXTRACT 0_2_0_result

### datファイルの作成

数値解析で剛体の運動表現したいことがよくある．
３次元で剛体の運動は，６自由度の運動で表現される．

$`t, x, y, z, \theta_x, \theta_y, \theta_z`$の順に並んでいる．

```data
0.0, 0.6666666421574411, -8.679523690686062e-05, 0., 0., 0., 0.00022580754704734906
0.01, 0.6666666360892092, -9.436993153716977e-05, 0., 0., 0., 0.0002085691971914809
0.02, 0.6666666288086627, -0.00010203474300463844, 0., 0., 0., 0.00018136038406231864
0.03, 0.666666620162383, -0.00010965062939973597, 0., 0., 0., 0.00014179367130820136
0.04, 0.6666666100078016, -0.00011703943260596665, 0., 0., 0., 8.706011682802684e-05
0.05, 0.6666665982273601, -0.0001239762852648563, 0., 0., 0., 1.3879082155409919e-05
```

'''

if True:
    time_data = []
    bodyA_data = {'x': [], 'y': [], 'angle': []}
    bodyB_data = {'x': [], 'y': [], 'angle': []}
    bodyC_data = {'x': [], 'y': [], 'angle': []}

    for i in range(0, 500):
        t = i * 0.01
        time_data.append(t)
        angles = robot.getAngles(t)
        c = (math.tanh((t-.5)*2.*math.pi)+1.)/2.
        for i in range(0, len(angles)):
            angles[i] *= c
        pos = robot.anglesToX(angles)
        bodyA_data['x'].append(pos[0][0])
        bodyA_data['y'].append(pos[0][1])
        bodyA_data['angle'].append(angles[0])
        bodyB_data['x'].append(pos[1][0])
        bodyB_data['y'].append(pos[1][1])
        bodyB_data['angle'].append(angles[1])
        bodyC_data['x'].append(pos[2][0])
        bodyC_data['y'].append(pos[2][1])
        bodyC_data['angle'].append(angles[2])

    def write_to_dat(filename, time_data, body_data):
        with open(filename, 'w') as f:
            f.write("# t, x, y, z, qx, qy, qz\n")  # Header
            for t, x, y, qz in zip(time_data, body_data['x'], body_data['y'], body_data['angle']):
                f.write(f"{t}, {x}, {y}, 0., 0., 0., {qz}\n")

    # Write data to .dat files
    write_to_dat("bodyA.dat", time_data, bodyA_data)
    write_to_dat("bodyB.dat", time_data, bodyB_data)
    write_to_dat("bodyC.dat", time_data, bodyC_data)

# --------------------------- datファイルを読み込んでグラフを表示する -------------------------- #

'''DOC_EXTRACT 0_2_0_result

### datファイルを読み込んでグラフを表示する

正しくdatファイルが作成できているか確認するために，datファイルを読み込んでグラフを表示する．

<img src="sample.png"  width="50%" height="50%">

'''


# Initialize lists to store angle and time data
if True:
    time_data = []
    x_data = []
    y_data = []
    qz_data = []

    fig = plt.figure()
    ax1 = fig.add_subplot(2, 1, 1)
    ax2 = fig.add_subplot(2, 1, 2)
    i = 0
    for name in ['bodyA', 'bodyB', 'bodyC']:
        with open(f'{name}.dat', 'r') as f:
            time_data.append([])
            x_data.append([])
            y_data.append([])
            qz_data.append([])
            for line in f.readlines():
                if line.startswith("#"):
                    continue  # Skip the line if it starts with #
                row = [float(x) for x in line.split(",")]
                time_data[i].append(row[0])
                x_data[i].append(row[1])
                y_data[i].append(row[2])
                qz_data[i].append(row[6])

        # Plot the datacode
        ax1.plot(time_data[i], qz_data[i], label='Angle'+str(i))
        ax2.plot(x_data[i], y_data[i], label='x,y'+str(i))
        i = i + 1

    # Label the axes
    ax1.set_ylabel('Angles')
    ax1.set_xlabel('Time')
    ax2.set_ylabel('y')
    ax2.set_xlabel('x')

    # Show the legend
    ax1.legend()

    # initial position of nodes
    for i in range(0, len(x_data)):
        print(i, ", x,y = ", x_data[i][0], " ", y_data[i][0])

    # save figure
    plt.savefig("sample.png")

    # Display the plot
    plt.show()
