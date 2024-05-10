import pyvista as pv
import numpy as np
import plotly.graph_objects as go

# Asegúrate de que este script se ejecute en un entorno que soporte la visualización interactiva, como Jupyter Notebook.

# Carga tu archivo .vtk aquí
filename = 'plane.vtk'  # Reemplaza 'tu_archivo.vtk' con la ruta de tu archivo.

# Leer el archivo .vtk
mesh = pv.read(filename)

# Verificar que el mesh es un UnstructuredGrid
if not isinstance(mesh, pv.UnstructuredGrid):
    raise TypeError("El archivo proporcionado no es un UnstructuredGrid.")

# Extracción de triángulos del UnstructuredGrid
cells = mesh.cells.reshape(-1, 4)  # Asumiendo que todos los triángulos vienen en bloques de 4 (1 valor de tamaño + 3 índices de vértices)
triangles = cells[:, 1:4]  # Ignoramos el primer valor de cada bloque, que es simplemente el tamaño de la célula (3 para triángulos)

# Si el script anteriormente no finalizaba, asegúrate de que 'mesh.cells' y 'mesh.points' tienen los datos esperados.

# Crear la figura Plotly para visualizar los triángulos
fig = go.Figure(data=[go.Mesh3d(x=mesh.points[:, 0],
                                y=mesh.points[:, 1],
                                z=mesh.points[:, 2],
                                i=triangles[:, 0],
                                j=triangles[:, 1],
                                k=triangles[:, 2],
                                opacity=0.5,
                                color='blue')])

# Configuración de la apariencia de la figura
fig.update_layout(scene=dict(xaxis_title='X',
                             yaxis_title='Y',
                             zaxis_title='Z',
                             xaxis=dict(visible=True),
                             yaxis=dict(visible=True),
                             zaxis=dict(visible=True)),
                  width=800, height=600)

# Mostrar la figura
fig.show()
