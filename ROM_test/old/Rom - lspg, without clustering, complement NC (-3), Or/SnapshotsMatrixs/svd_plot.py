import numpy as np
import matplotlib.pyplot as plt

# Importing the Kratos Library
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition

dir = ""
name_list = ["Phi_G_0.001"]

for name in name_list:

    SnapshotsMatrix = np.load(dir + name + ".npy")  

    _,s,_,_= RandomizedSingularValueDecomposition().Calculate(SnapshotsMatrix, 1e-12)

    valores_singulares = np.sort(s)[::-1]
    suma_acumulada = np.cumsum(valores_singulares)
    suma_acumulada_normalizada = suma_acumulada / suma_acumulada[-1]
    decaimiento = 1 - suma_acumulada_normalizada

    plt.figure(figsize=(10, 6))
    plt.semilogy(range(len(valores_singulares)), decaimiento, marker='x', markersize=1, linestyle='', color='b')
    plt.title('Singular Values Decay')
    plt.xlabel('Id singular Value')
    plt.ylabel('Decay')
    plt.grid(True, which="both", ls="--")
    plt.savefig(dir + name + '_Singular_Values_Decay.png', bbox_inches='tight')
    plt.show()
    plt.close()