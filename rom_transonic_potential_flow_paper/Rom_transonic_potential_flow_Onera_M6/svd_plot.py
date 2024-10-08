import numpy as np
import matplotlib.pyplot as plt

s = np.load('rom_data/SingularValuesVector.npy')

valores_singulares = np.sort(s)[::-1]

suma_acumulada = np.cumsum(valores_singulares)

suma_acumulada_normalizada = suma_acumulada / suma_acumulada[-1]

decaimiento = 1 - suma_acumulada_normalizada

plt.figure(figsize=(10, 6))
plt.semilogy(range(len(valores_singulares)), decaimiento, marker='x', markersize=5, linestyle='', color='b')
plt.title('Singular values decay')
plt.xlabel('Singular value index')
plt.ylabel('Decay')
plt.grid(True, which="both", ls="--")
plt.savefig('Singular_values_decay.png', bbox_inches='tight')
# plt.show()
plt.close()