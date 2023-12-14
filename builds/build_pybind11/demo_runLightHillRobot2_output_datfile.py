'''DOC_EXTRACT 0_2_0_result

### モーターの節の位置と角度の時間変化をdatファイルに出力する

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

# -------------------------------- datファイルの作成 -------------------------------- #

time_data = []
bodyA_data = {'x': [], 'y': [], 'angle': []}
bodyB_data = {'x': [], 'y': [], 'angle': []}
bodyC_data = {'x': [], 'y': [], 'angle': []}

num_of_data = 100
dt = 3*period / num_of_data

for i in range(0, num_of_data):
    t = i * dt
    time_data.append(t)
    angles = robot.getAngles(t)
    # c = (math.tanh((t-.5)*2.*math.pi)+1.)/2.
    # for i in range(0, len(angles)):
    #     angles[i] *= c
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
