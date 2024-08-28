import matplotlib.pyplot as plt
import numpy as np
import random

mu_values = np.load('mu.npy')
mu_values =  [mu.tolist() for mu in mu_values]

data_name = ["1.025, 0.72.dat", "1.0, 0.75.dat", "1.0, 0.73.dat"]

mach_train = [mu[1] for mu in mu_values]
alpha_train = [mu[0] for mu in mu_values]

# Ejemplo de datos (deberás reemplazar esto con tus datos reales)
np.random.seed(0)
# mach_train = np.random.rand(10)
# alpha_train = np.random.rand(10)
mach_test = np.random.rand(10)
alpha_test = np.random.rand(10)

# Crear figura y ejes
fig = plt.figure(figsize=(12, 12))

# Posición y tamaño del gráfico central (ajustable)
central_plot_position = [0.3, 0.3, 0.4, 0.4]  # [left, bottom, width, height]
ax = fig.add_axes(central_plot_position)
ax.scatter(mach_train, alpha_train, color='blue', label='Train')
ax.scatter(mach_test, alpha_test, color='red', label='Test')
ax.set_xlabel('Mach')
ax.set_ylabel('Alpha')
ax.legend()

# Posiciones relativas de los subplots
positions = [
    [0.08, 0.75, 0.2, 0.2],
    [0.72, 0.75, 0.2, 0.2],
    [0.08, 0.05, 0.2, 0.2],
    [0.72, 0.05, 0.2, 0.2]
]

# Dibujar subplots
for i, filename in enumerate(data_name):
    pos = positions[i]
    sub_ax = fig.add_axes(pos)
    x  = np.loadtxt(filename, usecols=(0,))
    cp = np.loadtxt(filename, usecols=(3,))
    sub_ax.plot( x, -cp, '.')
    sub_ax.set_title(f'Train {filename}')
    sub_ax.set_xlabel('x')
    sub_ax.set_ylabel('Cp')
    
    # Conectar con líneas o flechas
    arrowprops = dict(arrowstyle="->", color='black', lw=1)
    ax.annotate('',
                xy=(mach_train[i], alpha_train[i]), 
                xycoords='data',
                xytext=(pos[0] + pos[2]/2, pos[1] + pos[3]/2), 
                textcoords='figure fraction',
                arrowprops=arrowprops)

# Mostrar el gráfico
plt.show()
