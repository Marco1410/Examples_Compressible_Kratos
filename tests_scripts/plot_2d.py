import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

mu_train = np.load('mu_train.npy')
mu_train =  [mu.tolist() for mu in mu_train]

mu = mu_train[5]

# Cargar el archivo .vtk
mesh = pv.read(f"Results/vtk_output_FOM_Fit5/MainModelPart_Wing_0_1.vtk")

# Extraer la geometría superficial
surface_mesh = mesh.extract_geometry()

# Extraer los puntos y las celdas
points = surface_mesh.points
cells = surface_mesh.faces.reshape(-1, 4)[:, 1:4]

# Extraer los valores de la variable 'PRESSURE_COEFFICIENT'
pressure_coeff = surface_mesh.point_data['PRESSURE_COEFFICIENT']

# Filtrar puntos y celdas donde y > 0
mask = points[:, 1] > 0
filtered_points = points[mask]
filtered_pressure_coeff = pressure_coeff[mask]

# Crear un mapeo de los índices originales a los nuevos índices después del filtrado
index_map = {old_idx: new_idx for new_idx, old_idx in enumerate(np.where(mask)[0])}

# Reindexar las celdas para que correspondan a los índices en filtered_points
filtered_cells = []
for cell in cells:
    if all(mask[cell]):  # Solo conservar celdas donde todos los puntos tienen y > 0
        filtered_cells.append([index_map[idx] for idx in cell])

# Crear el gráfico 2D con interpolación de colores
plt.figure(figsize=(10, 8))
triang = tri.Triangulation(filtered_points[:, 0], filtered_points[:, 1], filtered_cells)
plt.tripcolor(triang, filtered_pressure_coeff, cmap='viridis')

# Configurar el gráfico
plt.colorbar(label='PRESSURE_COEFFICIENT')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Distribución de PRESSURE_COEFFICIENT para Y > 0')
plt.grid(True)
plt.axis('equal')

# Mostrar el gráfico
plt.show()
