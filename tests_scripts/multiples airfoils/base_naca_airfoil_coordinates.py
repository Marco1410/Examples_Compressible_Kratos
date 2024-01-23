import numpy as np

def coordinates_new_shape(x,y,z):
# Definir el perfil base: Naca Serie 4 digitos
# Longitud de la cuerda c = 1
# n Número de puntos en el perfil
# x vector de posiciones a lo largo de la cuerda
# Definir el rango de valores para el espesor máximo y curvatura
# naca xxxx posición del espesor máximo = 30%
# Ejemplo
# tmax_range = [0.10, 0.15] # Espesor máximo 
# m_range    = [0.0, 0.03] # Máxima curvatura
# p_range    = [0.25, 0.35] # Posicion de máxima curvatura

    filename = "airfoil_parameters.dat"
    tmax = np.loadtxt(filename,usecols=(0,))
    m = np.loadtxt(filename,usecols=(1,))
    p = np.loadtxt(filename,usecols=(2,))

    # Calcular linea de curvatura actual
    if x <= p:
        yc = (m/p**2)*(2*p*x-x**2)
        dyc = (2*m/p**2)*(p-x)
    else: 
        yc = (m/(1-p)**2)*((1-2*p)+2*p*x-x**2)
        dyc = (2*m/(1-p)**2)*(p-x)
    # Calcular el espesor máximo para el perfil actual
    yt = 5 * tmax * (0.2969*np.sqrt(x) - 0.1260*(x) - 0.3516*(x)**2 + 0.2843*(x)**3 - 0.1036*(x)**4)
    # Calcular las coordenadas del perfil actual
    if y > 0:
        xp = x - yt * np.sin(np.arctan(dyc)) 
        yp = yc + yt * np.cos(np.arctan(dyc))
    else:
        xp = x + yt * np.sin(np.arctan(dyc))
        yp = yc - yt * np.cos(np.arctan(dyc))
    return xp, yp, z
