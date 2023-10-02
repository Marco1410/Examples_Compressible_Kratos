import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
  
# Generate random 3D data points
x = np.loadtxt("parameters_list.dat",usecols=(0,))
y = np.loadtxt("parameters_list.dat",usecols=(1,))
z = np.loadtxt("parameters_list.dat",usecols=(2,))
data = np.array([x, y, z]).T
  
# Define mathematical function for curve fitting
def func(xy, a, b, c, d, e):
    x, y = xy
    return a /(1 + np.exp(-(b*x+c*y+d))) + e

# Perform curve fitting
popt, pcov = curve_fit(func, (x, y), z)
  
# Print optimized parameters
print(popt)
  
# Create 3D plot of the data points and the fitted curve
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, z, color='blue')
x_range = np.linspace(np.amin(x), np.amax(x), 50)
y_range = np.linspace(np.amin(y), np.amax(y), 50)

X, Y = np.meshgrid(x_range, y_range)
Z = func((X, Y), *popt)

ax.plot_surface(X, Y, Z, color='red', alpha=0.5)
ax.set_xlabel('Alpha')
ax.set_ylabel('Mach')
ax.set_zlabel('Upwind')
plt.show()