import numpy as np
import matplotlib.pyplot as plt

# Definir los parámetros de la función
A, B, C, D = 1.2, 0.4, 1.2, 10

x_shock = 0.6

# Definir la función combinada para Cp
def cp_function(x):
    if x < x_shock:
        return A * np.sin(np.pi * x)
    else:
        return B + C * np.tanh(D * (x - x_shock))

# Crear un rango de valores de x
x_values = np.linspace(0, 1, 100)

# Calcular los valores de Cp
cp_values = np.array([cp_function(x) for x in x_values])

# Graficar la función
plt.plot(x_values, cp_values)
plt.title('Modelo del Coeficiente de Presión (Cp) en un Flujo Transónico')
plt.xlabel('Posición a lo largo del perfil (x)')
plt.ylabel('Coeficiente de Presión (Cp)')
plt.grid(True)
plt.show()
