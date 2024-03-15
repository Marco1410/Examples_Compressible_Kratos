import numpy as np
import matplotlib.pyplot as plt

# Asumiendo que M es tu matriz original.
# Por ejemplo, M podría ser una matriz aleatoria para fines de demostración.
# Reemplaza esto con tu propia matriz.
M = np.load("M_full.npy")
n, m = M.shape  # n filas, m columnas

# Paso 1: Calcular SVD
U  = np.load("u.npy")  
s  = np.load("s.npy")
V = np.load("v.npy")
# vt = Vt.T


# Paso 2 y 3: Reconstruir M y calcular el error usando de 1 a m modos
errores = []
for i in range(1, m+1):
    S_i = np.diag(s[:i])  # Crear la matriz diagonal de valores singulares para i modos
    aux = (S_i @ Vt[:i, :])
    M_reconstruida = U[:, :i] @ aux  # Reconstruir M usando i modos
    error = np.linalg.norm(M - M_reconstruida, 'fro')  # Calcular el error de reconstrucción (norma Frobenius)
    errores.append(error)

# Paso 4: Graficar el error de reconstrucción
plt.figure(figsize=(10, 6))
plt.plot(range(1, m+1), errores, marker='o', linestyle='', color='b')
plt.title('Error de Reconstrucción de M usando de 1 a m Modos')
plt.xlabel('Número de Modos')
plt.ylabel('Error de Reconstrucción (Norma Frobenius)')
plt.grid(True)
plt.show()
