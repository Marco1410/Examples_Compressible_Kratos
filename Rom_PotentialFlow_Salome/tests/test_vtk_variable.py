import pyvista as pv
import plotly.graph_objects as go

filename = 'test.vtk'  # Reemplaza 'tu_archivo.vtk' con la ruta a tu archivo.

mesh = pv.read(filename)

# Asumimos que tienes al menos una variable escalar o vectorial asociada a los puntos o células.
# Por ejemplo, una variable escalar llamada 'PRESSURE_COEFFICIENT' y/o un campo vectorial llamado 'velocidad'.

# Para visualizar una variable escalar (ej. 'PRESSURE_COEFFICIENT') como colores en la malla:
if 'PRESSURE_COEFFICIENT' in mesh.point_data:
    PRESSURE_COEFFICIENT = mesh.point_data['PRESSURE_COEFFICIENT']
    escalar = PRESSURE_COEFFICIENT
elif 'PRESSURE_COEFFICIENT' in mesh.cell_data:
    PRESSURE_COEFFICIENT = mesh.cell_data['PRESSURE_COEFFICIENT']
    escalar = mesh.point_data_to_cell_data().cell_data['PRESSURE_COEFFICIENT']
else:
    escalar = None

# Preparar la visualización
cells = mesh.cells.reshape(-1, 4)  # Asumiendo que todos los triángulos vienen en bloques de 4 (1 valor de tamaño + 3 índices de vértices)
triangles = cells[:, 1:4]

fig = go.Figure(data=[go.Mesh3d(x=mesh.points[:, 0],
                                y=mesh.points[:, 1],
                                z=mesh.points[:, 2],
                                i=triangles[:, 0],
                                j=triangles[:, 1],
                                k=triangles[:, 2],
                                intensity=escalar,  # Usa la variable escalar aquí
                                colorscale='Viridis',  # Puedes elegir cualquier mapa de colores disponible
                                opacity=0.5)])

# Opciones adicionales para mejorar la visualización
fig.update_layout(scene=dict(xaxis_title='X',
                             yaxis_title='Y',
                             zaxis_title='Z',
                             xaxis=dict(visible=True),
                             yaxis=dict(visible=True),
                             zaxis=dict(visible=True)),
                  width=1024, height=780)

# Para visualizar la malla, hemos añadido la geometría directamente.
# Las capacidades de visualización de streamlines en Plotly son limitadas.
# Para un análisis detallado de campos vectoriales, considera usar software específico como ParaView.

fig.show()

# Guardar la figura como una imagen PNG
# fig.write_image("test.png")

