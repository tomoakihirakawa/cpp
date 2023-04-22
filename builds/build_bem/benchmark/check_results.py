import json
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

# Load JSON data
json_file_path = os.path.expanduser("~/BEM/Kramer2021_0d3/result.json")
with open(json_file_path, 'r') as file:
    data = json.load(file)

# Extract the data
float_COM = np.array(data['float_COM'])
float_COM = (float_COM - 0.9) / 0.074*2.5
time = np.array(data['time'])
z = float_COM[:, 2]

# Load text data
txt_file_path = './Datafile/Experimental_results/01D_Measured1_Normalized.txt'

# txt_data = pd.read_csv(txt_file_path, sep='\s+', engine='python')
txt_data = pd.read_csv(txt_file_path, sep='\s+',
                       engine='python', header=None, skiprows=1)
txt_data.columns = ['t/Te0 [-]', 'x3/H_{0,m} [-]',
                    'WG1/H_{0,m} [-]', 'WG2/H_{0,m} [-]', 'WG3/H_{0,m} [-]']

print(txt_data.columns)

# Extract text data
txt_time = txt_data['t/Te0 [-]']
txt_z = txt_data['x3/H_{0,m} [-]']

# Plot JSON data
plt.figure()
plt.plot(time, z, label="JSON data")

# Plot text data
plt.plot(txt_time, txt_z, label="Text data")

plt.xlabel("Time")
plt.ylabel("Value")
plt.title("Float Position Data")
plt.legend()
plt.show()
