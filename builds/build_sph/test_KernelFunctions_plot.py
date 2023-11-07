import matplotlib.pyplot as plt

# Read the data from the text file
x_values = []
bspline3_values = []
bspline4_value = []
bspline5_values = []
wendland_values = []
spike_values = []

with open("data.txt", "r") as file:
    for line in file:
        x, bspline3, bspline4, bspline5, wendland, spike = map(float, line.split())
        x_values.append(x)
        bspline3_values.append(bspline3)
        bspline4_value.append(bspline4)
        bspline5_values.append(bspline5)
        wendland_values.append(wendland)
        spike_values.append(spike)

# Plotting the data
plt.figure(figsize=(10, 6))

plt.plot(x_values, bspline3_values, label="B-spline 3", linestyle="-", marker="o")
plt.plot(x_values, bspline4_value, label="B-spline 4", linestyle="-.", marker="^")
plt.plot(x_values, bspline5_values, label="B-spline 5", linestyle="--", marker="x")
plt.plot(x_values, wendland_values, label="Wendland", linestyle=":", marker="s")
plt.plot(x_values, spike_values, label="Spike", linestyle="-.", marker="d")

plt.xlabel("x")
plt.ylabel("Kernel Value")
plt.title("Kernel Functions")
plt.legend()

plt.grid(True)
plt.show()
