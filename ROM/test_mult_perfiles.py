import os
from scipy.stats import qmc
import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt

def plot_mu_values(mu,name):
    # Creamos la figura
    fig = plt.figure()
    # Creamos el plano 3D
    ax1 = fig.add_subplot(111, projection='3d')
    # Definimos los datos
    mu_t = np.zeros(len(mu))
    mu_m = np.zeros(len(mu))
    mu_p = np.zeros(len(mu))
    for i in range(len(mu)):
        mu_t[i] = mu[i][0]
        mu_m[i] = mu[i][1]
        mu_p[i] = mu[i][2]
    # Agregamos los puntos en el plano 3D
    ax1.scatter(mu_t, mu_m, mu_p, c='r', marker='o')
    ax1.set_xlabel('t')
    ax1.set_ylabel('m')
    ax1.set_zlabel('p')    
    plt.title('Foil Parameters')
    plt.grid()
    plt.savefig(name)

def get_multiple_params_by_Halton(number_of_values,tmax,m,p):
    sampler = qmc.Halton(d=3)
    sample = sampler.random(number_of_values)
    mu = []
    values = qmc.scale(sample,[tmax[0],m[0],p[0]],[tmax[1],m[1],p[1]])
    for i in range(number_of_values):
        mu.append([values[i,0], values[i,1], values[i,2]])
    return mu

def new_foil(x,tmax,m,p):
    xp = np.array(np.zeros(2*len(x)))
    yp = np.array(np.zeros(2*len(x)))
    for i in range(len(x)):
        # Calcular linea de curvatura actual
        if x[i] <= p:
            yc = (m/p**2)*(2*p*x[i]-x[i]**2)
            dyc = (2*m/p**2)*(p-x[i])
        else: 
            yc = (m/(1-p)**2)*((1-2*p)+2*p*x[i]-x[i]**2)
            dyc = (2*m/(1-p)**2)*(p-x[i])
        # Calcular el espesor máximo para el perfil actual
        yt = 5 * tmax * (0.2969*np.sqrt(x[i]) - 0.1260*(x[i]) - 0.3516*(x[i])**2 + 0.2843*(x[i])**3 - 0.1036*(x[i])**4)
        # Calcular las coordenadas del perfil actual
        xp[i] = x[i] - yt * np.sin(np.arctan(dyc)) 
        yp[i] = yc + yt * np.cos(np.arctan(dyc))
        xp[(2*len(x)-1)-i] = x[i] + yt * np.sin(np.arctan(dyc))
        yp[(2*len(x)-1)-i] = yc - yt * np.cos(np.arctan(dyc))
    return xp, yp

# Definir el perfil base: Naca Serie 4 digitos
# Longitud de la cuerda c = 1
n = 150    # Número de puntos en el perfil
x = np.linspace(0, 1, n)    # Posiciones a lo largo de la cuerda

# Plot base airfoil ...............
fig, ax = plt.subplots()
fig.set_figwidth(15.0)
fig.set_figheight(6.8)

# Definir el rango de valores para el espesor máximo y curvatura
# naca xxxx posición del espesor máximo = 30%
tmax_range = [0.10, 0.15] # Espesor máximo 
m_range    = [0.0, 0.03] # Máxima curvatura
p_range    = [0.25, 0.35] # Posicion de máxima curvatura

mu_values = get_multiple_params_by_Halton(15,tmax_range,m_range,p_range)
plot_mu_values(mu_values,"parameters.png")

# Calcular el perfil modificado para cada conjunto de datos
for mu in mu_values:
    tmax,m,p = mu
    xp,yp = new_foil(x,tmax,m,p)
    ax.plot( xp, yp, "+", markersize=3)

ax.grid()
plt.axis('auto')
plt.xlabel("x")
plt.ylabel("y")
plt.tight_layout()
fig.savefig(os.path.join("airfoils.png"))