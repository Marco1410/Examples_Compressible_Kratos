import dask.array as da
import numpy as np
import gc
from dask.distributed import Client, LocalCluster



if __name__ == '__main__':
    # Configurar el clúster local
    cluster = LocalCluster(n_workers=8, threads_per_worker=1, memory_limit='2GB')
    client = Client(cluster)

    mu_values = np.arange(0, 300, 1, dtype=int)
    chunks = []

    for mu in mu_values:

        # Convertir a un arreglo Dask
        dask_array = da.from_array(np.load("data.npy"), chunks=(1000, 10))
        chunks.append(dask_array)

        # Liberar la memoria del arreglo de NumPy
        gc.collect()

    # Concatenar todos los arreglos de Dask
    ResidualsSnapshotsMatrix = da.concatenate(chunks, axis=1)

    # Persistir en el clúster para reducir el uso de memoria local
    ResidualsSnapshotsMatrix = ResidualsSnapshotsMatrix.persist()

    # Liberar la lista de fragmentos
    del chunks
    gc.collect()

    # Realiza operaciones adicionales
    print(ResidualsSnapshotsMatrix.shape)

    # Liberar memoria
    del ResidualsSnapshotsMatrix
    gc.collect()

    # Cerrar el cliente y el clúster cuando hayas terminado
    client.close()
    cluster.close()
