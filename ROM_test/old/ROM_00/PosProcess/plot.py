import numpy as np
import matplotlib.tri as tri
import matplotlib.pyplot as plt

snapshot = np.load(f"2.147843426588291, 0.7205444154266182.npy")
x = (np.loadtxt("nodes.dat", usecols=(1,))).flatten()
y = (np.loadtxt("nodes.dat", usecols=(2,))).flatten()

fig1, ax1 = plt.subplots()
ax1.set_aspect('equal')
tpc = ax1.tripcolor( x, y, snapshot[len(x)::,0], edgecolors='k')
fig1.colorbar(tpc)
ax1.set_title('tripcolor of Delaunay triangulation, flat shading')
plt.show()