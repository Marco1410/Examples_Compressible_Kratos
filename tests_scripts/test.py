import matplotlib.pyplot as plt
import numpy as np


mu_values = np.load('mu.npy')
mu_values =  [mu.tolist() for mu in mu_values]

data_name = ["1.025, 0.72.dat", "1.0, 0.75.dat", "1.0, 0.73.dat"]

# Ejemplo de datos (deberás reemplazar esto con tus datos reales)
mach = [mu[1] for mu in mu_values]
alpha = [mu[0] for mu in mu_values]

# Número de conjuntos de datos
num_datasets = len(mach)

# Crear figura y ejes
fig = plt.figure(figsize=(10, 10))

# Posiciones relativas de los subplots
positions = [
    (0.1, -1.6, 1, 1),
    (1.6, -1.6, 1, 1),
    (1.6, -0.1, 1, 1)
]

center_ax = fig.add_axes([0.3, 0.3, 0.2, 0.2])
# Dibujar subplots
for i in range(len(positions)):
    ax = center_ax.inset_axes(positions[i])
    x  = np.loadtxt(data_name[i], usecols=(0,))
    cp = np.loadtxt(data_name[i], usecols=(3,))
    ax.plot(x,-cp,'.')
    ax.set_xlim(1, 1.01)
    ax.set_xticks([1, 1.01])
    ax.set_xticklabels(['-1.0', '1.5oo'])
    ax.set_ylim(1.0, 1.01)
    ax.set_yticks([1.0, 1.01])
    ax.set_yticklabels(['-1.0', '1.5oo'])
    # ax.set_xlabel('X')
    # ax.set_ylabel('Cp')
    ax.set_title(f'Mach={mach[i]}, Alpha={alpha[i]}')
    center_ax.indicate_inset_zoom(ax, edgecolor="grey")

# Dibujar gráfico central
center_ax.plot(mach, alpha, 'o')
center_ax.set_xlabel('Mach')
center_ax.set_ylabel('Alpha')
center_ax.set_title('Parameter Space')

# Mostrar el gráfico
plt.show()
