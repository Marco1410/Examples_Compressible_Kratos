import matplotlib.pyplot as plt
import numpy as np
import random

mu_train = np.load('mu_train.npy')
mu_train =  [mu.tolist() for mu in mu_train]
mach_train = [mu[1] for mu in mu_train]
alpha_train = [mu[0] for mu in mu_train]

mu_test = np.load('mu_test.npy')
mu_test =  [mu.tolist() for mu in mu_test]
mach_test = [mu[1] for mu in mu_test]
alpha_test = [mu[0] for mu in mu_test]

mu_validation = np.load('mu_validation.npy')
mu_validation =  [mu.tolist() for mu in mu_validation]
mach_validation = [mu[1] for mu in mu_validation]
alpha_validation = [mu[0] for mu in mu_validation]

# Elegir pares aleatorios del set train
# indices = random.sample(range(len(mach_train)), 4)
# print(indices)
indices = [51, 99, 48, 6]

# Dibujar subplots
# Crear figura y ejes
fig = plt.figure(figsize=(10, 10))

# Posición y tamaño del gráfico central (ajustable)
central_plot_position = [0.3, 0.3, 0.4, 0.4]  # [left, bottom, width, height]
ax = fig.add_axes(central_plot_position)  

regions = [
    ((min(mach_train), min(alpha_train)), (max(mach_train), max(alpha_train)), 'orange')    
]
for bottom_left, top_right, color in regions:
    rect = plt.Rectangle(bottom_left, top_right[0] - bottom_left[0], top_right[1] - bottom_left[1], 
                        facecolor=color, edgecolor=color, alpha=0.25)
    ax.add_patch(rect)

ax.plot(mach_train     , alpha_train     , 'bs', label='Train'     , markersize=4)
ax.plot(mach_test      , alpha_test      , 'rx', label='Test'      , markersize=7)
ax.plot(mach_validation, alpha_validation, 'g*', label='Validation', markersize=10)
# ax.set_xlim(0.699,0.751)
ax.set_xlabel('Mach')
ax.set_ylabel('Alpha')
ax.legend(bbox_to_anchor=(1.01, 0.9, 1.0, 0.1), loc='upper left', borderaxespad=0.0)
ax.set_title('Parameter space', size=13, pad=8)

# Posiciones relativas de los subplots
positions = [
    [0.08, 0.75, 0.2, 0.2],
    [0.75, 0.75, 0.2, 0.2],
    [0.08, 0.05, 0.2, 0.2],
    [0.75, 0.05, 0.2, 0.2]
]

# Dibujar subplots
I=II=III=IV=False
for idx in indices:
    filename = f'FOM_Skin_Data/{alpha_train[idx]}, {mach_train[idx]}.dat'
    if (mach_train[idx] <= np.min(mach_train) + ((np.max(mach_train)-np.min(mach_train))/2)):
        if (alpha_train[idx] <= np.min(alpha_train) + ((np.max(alpha_train)-np.min(alpha_train))/2)):
            if not III:
                pos = positions[2]
                III = True
            else: 
                pos = positions[0]
        else:
            if not I:
                pos = positions[0]
                I = True
            else: 
                pos = positions[2]
    else:
        if (alpha_train[idx] <= np.min(alpha_train) + ((np.max(alpha_train)-np.min(alpha_train))/2)):
            if not IV:
                pos = positions[3]
                IV = True
            else: 
                pos = positions[1]
        else:
            if not II:
                pos = positions[1]
                II = True
            else: 
                pos = positions[3]
    sub_ax = fig.add_axes(pos)
    x  = np.loadtxt(filename, usecols=(0,))
    cp = np.loadtxt(filename, usecols=(3,))
    sub_ax.plot( x, -cp, '.', markersize = 1.0)
    sub_ax.set_title(f'Angle={np.round(alpha_train[idx],3)}, Mach={np.round(mach_train[idx],3)}')
    sub_ax.set_ylim(-1.0,1.5)
    sub_ax.set_xlabel('x')
    sub_ax.set_ylabel('-Cp')
    
    # Conectar con líneas o flechas
    arrowprops = dict(arrowstyle="->", color='grey', lw=1)
    ax.annotate('',
                xy=(mach_train[idx], alpha_train[idx]), 
                xycoords='data',
                xytext=(pos[0] + pos[2]/2, pos[1] + pos[3]/2), 
                textcoords='figure fraction',
                arrowprops=arrowprops)

# Mostrar el gráfico
# plt.show()
plt.savefig('Mu_parameters_cp.pdf', dpi=400)
plt.close('all')
