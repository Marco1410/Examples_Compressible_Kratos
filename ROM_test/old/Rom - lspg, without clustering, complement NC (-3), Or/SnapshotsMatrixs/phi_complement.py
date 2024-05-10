from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

SnapshotsMatrix     = np.load("SnapshotsMatrix_nc_0.npy")
aux_SnapshotsMatrix = np.load("SnapshotsMatrix_conv_0.npy")
tolerances = [1e-3, 1e-4, 1e-5]

plt.figure(figsize=(10, 6))
colors = cm.viridis(np.linspace(0, 1, len(tolerances)))

phi_c,_,_,_= RandomizedSingularValueDecomposition().Calculate(aux_SnapshotsMatrix, 1e-12)

S_nc = SnapshotsMatrix
for k in range(5):
    S_nc = S_nc - (phi_c @ (phi_c.T @ S_nc))

for id, tolerance in enumerate(tolerances):

    phi_nc,_,_,_= RandomizedSingularValueDecomposition().Calculate(S_nc, tolerance)
    phi_G = np.concatenate((phi_c, phi_nc), axis=1)
    np.save("Phi_G_"+str(tolerance), phi_G)

    u,s,_,_= RandomizedSingularValueDecomposition().Calculate(phi_G, 1e-12)
    valores_singulares = np.sort(s)[::-1]
    suma_acumulada = np.cumsum(valores_singulares)
    suma_acumulada_normalizada = suma_acumulada / suma_acumulada[-1]
    decaimiento = 1 - suma_acumulada_normalizada
    plt.semilogy(range(len(valores_singulares)), decaimiento, label='tol: '+str(tolerance)+' modes: '+str(u.shape[1]), markersize=1, linestyle='', color=colors[id])

plt.title('Singular Values Decay')
plt.xlabel('Id singular Value')
plt.ylabel('Decay')
plt.grid(True, which="both", ls="--")
plt.savefig('Singular_Values_Decay.png', bbox_inches='tight')
plt.show()
plt.close()