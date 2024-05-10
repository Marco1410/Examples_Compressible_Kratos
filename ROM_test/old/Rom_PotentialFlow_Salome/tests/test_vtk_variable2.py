import pyvista as pv
import plotly.graph_objects as go

filename = 'plane.vtk'
mesh = pv.read(filename)

mesh.set_active_scalars('PRESSURE_COEFFICIENT')

if hasattr(mesh, 'cells'):
    cells = mesh.cells.reshape(-1, 4)
    is_triangle = cells[:, 0] == 3
    triangles = cells[is_triangle, 1:4]

    pressure_coefficient = mesh.point_data['PRESSURE_COEFFICIENT']
    mesh_fig = go.Mesh3d(x=mesh.points[:, 0], y=mesh.points[:, 1], z=mesh.points[:, 2],
                         i=triangles[:, 0], j=triangles[:, 1], k=triangles[:, 2],
                         intensity=pressure_coefficient,
                         colorscale='Jet', colorbar=dict(title='Pressure Coefficient'),
                         opacity=0.5)

    fig = go.Figure(data=[mesh_fig])

    # Configurar la vista inicial (opcional)
    camera = dict(
        up=dict(x=0, y=1, z=0),
        center=dict(x=0.003, y=0, z=0),
        eye=dict(x=0.003, y=0, z=0.025)  # Ajusta estos valores según sea necesario para tu visualización
    )

    # Botón para cambiar a la vista del plano XY
    fig.update_layout(
        updatemenus=[
            dict(
                type="buttons",
                direction="left",
                buttons=[
                    dict(
                        args=[{"scene.camera": camera}],
                        label="Vista Plano XY",
                        method="relayout"
                    )
                ],
                pad={"r": 10, "t": 10},
                showactive=True,
                x=0.11,
                xanchor="left",
                y=1.08,
                yanchor="top"
            ),
        ],
        scene=dict(
            xaxis_title='X', yaxis_title='Y', zaxis_title='Z',
            xaxis=dict(visible=True), yaxis=dict(visible=True),
            zaxis=dict(visible=True),
            camera=camera  # Establecer la cámara inicial
        ),
        width=1024,
        height=780,
        title="Naca 0012 PRESSURE_COEFFICIENT"
    )

    fig.show()
else:
    raise AttributeError("El objeto UnstructuredGrid no tiene el atributo 'cells'.")
