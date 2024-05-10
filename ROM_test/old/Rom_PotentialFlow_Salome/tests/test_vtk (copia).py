import pyvista as pv
from dash import Dash, dcc, html
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
                         colorscale='Viridis', colorbar=dict(title='Pressure Coefficient'),
                         opacity=1.0)

    fig = go.Figure(data=[mesh_fig])

# Update plot sizing
fig.update_layout(
    width=1024,
    height=780,
    autosize=False,
    margin=dict(t=0, b=0, l=0, r=0),
    scene_camera=dict(eye=dict(x=1.0, y=3.0, z=3.0)),
    template="plotly_dark",
)

# Update 3D scene options
fig.update_scenes(
    aspectratio=dict(x=4, y=2, z=0.5),
    aspectmode="data"
)
# fig.show()

app = Dash()
app.layout = html.Div([
    dcc.Graph(figure=fig)
])

app.run_server(debug=True, use_reloader=True)  # Turn off reloader if inside Jupyter
