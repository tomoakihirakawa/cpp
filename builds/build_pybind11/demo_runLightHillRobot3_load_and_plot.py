'''DOC_EXTRACT 0_3_0_result

### 作成したdatファイルを読み込んで確認する

![sample.png](sample.png)

'''

import math
import numpy as np
import matplotlib.pyplot as plt

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
