import pyvista as pv
import numpy as np
import plotly.graph_objects as go

filename = 'test.vtk'  # Reemplaza con la ruta a tu archivo
mesh = pv.read(filename)

mesh.set_active_scalars('PRESSURE_COEFFICIENT')

# Generar contornos
min_val, max_val = mesh.get_data_range()
num_contours = 10
contours = mesh.contour(isosurfaces=np.linspace(min_val, max_val, num_contours))

# Visualizar la geometría de la malla
geom_x, geom_y, geom_z = mesh.points.T
mesh_fig = go.Mesh3d(x=geom_x, y=geom_y, z=geom_z, opacity=0.1, color='gray')

# Preparar las listas de coordenadas para las líneas de contorno
x, y, z = [], [], []

for i in range(0, len(contours.lines), 3):
    point_ids = contours.lines[i+1:i+3]  # Los ids de los puntos están en las posiciones i+1 e i+2
    for pid in point_ids:
        x.append(contours.points[pid][0])
        y.append(contours.points[pid][1])
        z.append(contours.points[pid][2])
    # Añadir None después de cada segmento para no conectar las líneas entre sí
    x.append(None)
    y.append(None)
    z.append(None)

# Crear figura Plotly para visualizar la malla y las líneas de contorno
fig = go.Figure(data=[mesh_fig])

# Añadir las líneas de contorno
fig.add_trace(go.Scatter3d(x=x, y=y, z=z, mode='lines', line=dict(color='blue', width=2)))

# Configurar la apariencia de la figura
fig.update_layout(title="Geometría y Líneas de Contorno de PRESSURE_COEFFICIENT",
                  scene=dict(xaxis_title='X', yaxis_title='Y', zaxis_title='Z',
                             xaxis=dict(visible=True), yaxis=dict(visible=True), zaxis=dict(visible=True)),
                  width=800, height=600)

# Mostrar la figura
fig.show()
