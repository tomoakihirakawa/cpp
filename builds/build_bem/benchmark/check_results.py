import json
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

# Load JSON data
file = "~/BEM/Kramer2021_H00d03_small/result.json"
# file = "~/BEM/Kramer2021_H00d03/result.json"
json_file_path = os.path.expanduser(file)
with open(json_file_path, 'r') as file:
    data = json.load(file)

# Extract the data
float_COM = np.array(data['float_COM'])
Te0 = 0.76
time = np.array(data['time'])/Te0
H0 = 30/1000
z = (float_COM[:, 2] - 900 / 1000) / H0

# Load text data
txt_file_path = './Datafile/Experimental_results/01D_CI95_Normalized.txt'
exp_data = pd.read_csv(txt_file_path, sep='\s+',
                       engine='python', header=None, skiprows=1)
exp_data.columns = ['t/Te0 [-]',
                    'x3/H_{0,m} (mean) [-]',
                    'Lower 95% CI bound [-]',
                    'Upper 95% CI bound [-]']

print(exp_data.columns)

# Extract text data
txt_time = exp_data['t/Te0 [-]']
txt_z = exp_data['x3/H_{0,m} (mean) [-]']


# Plot JSON data
plt.figure()
plt.plot(time, z, 'go-', label="JSON data")


# Plot text data
plt.plot(txt_time, txt_z, label="Text data")

plt.xlabel("Time")
plt.ylabel("Value")

plt.ylim([-1.1, 1.1])
plt.xlim([0, 2])

plt.title("Float Position Data")
plt.legend()
plt.show()
