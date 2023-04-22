import numpy as np
import pandas as pd
import os
import json
import matplotlib.pyplot as plt

notebook_directory = "./"
data01DLPF4 = pd.read_csv(os.path.join(
    notebook_directory, "Datafile/Numerical_results/LPF4/01D_LPF4.txt"), sep="\s+", header=None)
data01DFNPF1 = pd.read_csv(os.path.join(
    notebook_directory, "Datafile/Numerical_results/FNPF1/01D_FNPF1.txt"), sep="\s+", header=None)
data01DMeasured1Raw = pd.read_csv(os.path.join(
    notebook_directory, "Datafile/Experimental_results/01D_CI95_Normalized.txt"), sep="\s+", header=None)
data03DMeasured1Raw = pd.read_csv(os.path.join(
    notebook_directory, "Datafile/Experimental_results/03D_Measured1_Normalized.txt"), sep="\s+", header=None)
data05DMeasured1Raw = pd.read_csv(os.path.join(
    notebook_directory, "Datafile/Experimental_results/05D_Measured1_Normalized.txt"), sep="\s+", header=None)

name = os.path.expanduser("~/BEM/Kramer2021_0d3/result.json")
with open(name) as f:
    data_json = json.load(f)

data = np.array(data_json["float_COM"])
time = np.array(data_json["time"])
Z = data[:, 2]
H0 = 30 / 1000

x1 = data01DMeasured1Raw.iloc[1:, 0] / 0.76
y1 = data01DMeasured1Raw.iloc[1:, 1] / H0

x2 = (time - 0.011) / 0.76
y2 = (Z - (-34.8 + 900) / 1000) / H0

plt.plot(x1, y1, 'o', color='darkgreen', alpha=0.6, label="Data 1")
plt.plot(x2, y2, color='red', label="Data 2")

plt.xlim(0, 2)
plt.ylim(-1, 1)
plt.xlabel("x-axis")
plt.ylabel("y-axis")
plt.title("Scientific Plot")
plt.grid(True)
plt.legend()

plt.gca().xaxis.set_ticks(np.arange(0, 4, 0.2))
plt.gca().yaxis.set_ticks(np.arange(-1, 1, 0.5))

plt.show()
