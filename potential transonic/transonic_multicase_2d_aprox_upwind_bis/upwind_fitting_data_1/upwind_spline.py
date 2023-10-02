import numpy as np
from scipy.interpolate import Rbf
import matplotlib.pyplot as plt

# Generate random 3D data points
x = np.loadtxt("parameters_list.dat",usecols=(0,))
y = np.loadtxt("parameters_list.dat",usecols=(1,))
z = np.loadtxt("parameters_list.dat",usecols=(2,))
x_range = np.linspace(np.amin(x), np.amax(x), 50)
y_range = np.linspace(np.amin(y), np.amax(y), 50)
X, Y = np.meshgrid(x_range, y_range)
# Fit a radial basis function model        
# 'multiquadric','inverse','gaussian','linear'
# 'cubic','quintic','thin_plate'
rbf = Rbf(x, y, z, function="multiquadric", epsilon=0.001)
Z = rbf(X, Y)

# Plot the original data and the fitted function
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.scatter(x, y, z, color='blue')
ax.plot_surface(X, Y, Z, color='red', alpha=0.5)
ax.set_xlabel('Alpha')
ax.set_ylabel('Mach')
ax.set_zlabel('Upwind')
plt.show()
