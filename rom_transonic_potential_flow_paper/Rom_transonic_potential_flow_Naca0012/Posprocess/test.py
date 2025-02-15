import numpy as np
import matplotlib.pyplot as plt

Cp_FOM_3D = np.loadtxt('FOM_0.45, 0.73.dat', usecols=(3,))
Cp_ROM_3D = np.loadtxt('ROM_0.45, 0.73.dat', usecols=(3,))

plt.figure(figsize=(6, 5))
plt.scatter(Cp_FOM_3D, Cp_ROM_3D, c='w', alpha=0.5, edgecolors='r', label="Data")
val_max = np.max(Cp_FOM_3D)
val_min = np.min(Cp_FOM_3D)
plt.plot([val_min, val_max], [val_min, val_max], 'k--', linewidth=1, label="y=x")
plt.xlabel("$C_p$ FOM")
plt.ylabel("$C_p$ ROM")
plt.title("Comparación ROM vs FOM")
plt.grid()
plt.legend()
plt.show()
