import matplotlib.pyplot as plt
import numpy as np

# Load data
himmelblau_data = np.loadtxt('Himmelblau.dat')
himmelblau_broyden_data = np.loadtxt('Himmelblau_Broyden.dat')
himmelblau_newton_data = np.loadtxt('Himmelblau_Newton.dat')  # load new data

# Create a 2D grid for the contour plot
x = np.unique(himmelblau_data[:,0])
y = np.unique(himmelblau_data[:,1])
X, Y = np.meshgrid(x, y)
Z = himmelblau_data[:,2].reshape(len(x), len(y))

# Contour plot
plt.figure(figsize=(8, 6))

# Define contour levels
levels = []
for i in range(-100, 1000, 5):
    levels.append(i)

contourf = plt.contourf(Y, X, Z, levels=levels, cmap='coolwarm')  # change colormap here
contour = plt.contour(Y, X, Z, levels=levels, colors='black', linewidths=0.1)  # make lines thinner

# Scatter plot
plt.plot(himmelblau_broyden_data[:,0], himmelblau_broyden_data[:,1], 'x-', color='orange', label='Broyden points')  # plot Broyden data with orange color and circular marker
plt.plot(himmelblau_newton_data[:,0], himmelblau_newton_data[:,1], 'x-', color='magenta', label='Newton points')  # plot Newton data with magenta color and circular marker

plt.title('Contour plot of Himmelblau function with Broyden and Newton Method optimization points')
plt.xlim(-6, 6)
plt.ylim(-6, 6)
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.colorbar(contourf, label='F(x, y)')  # add a colorbar

plt.savefig('output.png', dpi=300)
plt.show()
