import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
  
# Generate random 3D data points
x = np.random.random(100)
y = np.random.random(100)
z = np.sin(x * y) + np.random.normal(0, 0.1, size=100)
data = np.array([x, y, z]).T
  
# Define mathematical function for curve fitting
def func(xy, a, b, c, d, e, f):
    x, y = xy
    return a + b*x + c*y + d*x**2 + e*y**2 + f*x*y
  
# Perform curve fitting
popt, pcov = curve_fit(func, (x, y), z)
  
# Print optimized parameters
print(popt)
  
# Create 3D plot of the data points and the fitted curve
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, z, color='blue')
x_range = np.linspace(0, 1, 50)
y_range = np.linspace(0, 1, 50)
X, Y = np.meshgrid(x_range, y_range)
Z = func((X, Y), *popt)
ax.plot_surface(X, Y, Z, color='red', alpha=0.5)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()