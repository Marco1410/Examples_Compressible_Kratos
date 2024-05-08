import numpy as np
import matplotlib.pyplot as plt

# Importing the Kratos Library
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition

dir = "SnapshotsMatrixs/"
name_list = ["SnapshotsMatrix_conv_0",
             "SnapshotsMatrix_full_0"]

for name in name_list:

    SnapshotsMatrix = np.load(dir + name + ".npy")  

    _,s,_,_= RandomizedSingularValueDecomposition().Calculate(SnapshotsMatrix, 1e-12)

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
    plt.savefig(dir + name + '_Decaimiento_valores_singulares.png', bbox_inches='tight')
    # plt.show()
    plt.close()

    plt.figure(figsize=(10, 6))
    plt.plot(range(len(valores_singulares)), valores_singulares, marker='x', markersize=1, linestyle='', color='b')
    plt.title('Valores singulares')
    plt.xlabel('Índice del Valor Singular')
    plt.ylabel('Magnitud')
    plt.grid(True, which="both", ls="--")
    plt.savefig(dir + name + '_Valores_singulares.png', bbox_inches='tight')
    # plt.show()
    plt.close()