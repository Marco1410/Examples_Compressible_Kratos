import numpy as np
from scipy.interpolate import Rbf
import matplotlib.pyplot as plt

# Generate random 3D data points
x = np.linspace(-5, 5, 100)
y = np.linspace(-5, 5, 100)
X, Y = np.meshgrid(x, y)
Z = np.cos(np.sqrt(X**2 + Y**2))
# Fit a radial basis function model
rbf = Rbf(X, Y, Z, function="quintic")
Z_pred = rbf(X, Y)

# Plot the original data and the fitted function
fig = plt.figure()
ax = fig.add_subplot(projection="3d")
ax.plot_surface(X, Y, Z)
ax.plot_surface(X, Y, Z_pred)
plt.show()
