import numpy as np
import matplotlib.pyplot as plt

# Importing the Kratos Library
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition


SnapshotsMatrix_c  = np.load("SnapshotsMatrix_conv.npy")  
SnapshotsMatrix_nc = np.load("SnapshotsMatrix_not_conv.npy")  

# print(np.concatenate((SnapshotsMatrix_c,SnapshotsMatrix_nc),axis=1).shape)

phi_c,_,_,_= RandomizedSingularValueDecomposition().Calculate(SnapshotsMatrix_c, 1e-12)

S_nc = SnapshotsMatrix_nc
for k in range(5):
    S_nc = S_nc - (phi_c @ (phi_c.T @ S_nc))
phi_nc,_,_,_= RandomizedSingularValueDecomposition().Calculate(S_nc, 1e-6)

phi = np.linalg.svd(np.c_[phi_c, phi_nc], full_matrices=False)[0]

_,s,_,_= RandomizedSingularValueDecomposition().Calculate(phi, 1e-6)

# valores singulares ordenados de mayor a menor.
valores_singulares = np.sort(s)[::-1]

# Calcula la suma acumulada de los valores singulares
suma_acumulada = np.cumsum(valores_singulares)

# Normaliza dividiendo por la suma total de los valores singulares
suma_acumulada_normalizada = suma_acumulada / suma_acumulada[-1]

# Calcula 1 - suma acumulada normalizada para obtener la información no capturada o decaimiento
decaimiento = 1 - suma_acumulada_normalizada

plt.figure(figsize=(10, 6))
plt.semilogy(range(len(valores_singulares)), decaimiento, marker='x', markersize=1, linestyle='', color='b')
plt.title('Decaimiento valores singulares')
plt.xlabel('Índice del Valor Singular')
plt.ylabel('Decaimiento')
plt.grid(True, which="both", ls="--")
# plt.savefig('Decaimiento_valores_singulares.png', bbox_inches='tight')
plt.show()
plt.close()