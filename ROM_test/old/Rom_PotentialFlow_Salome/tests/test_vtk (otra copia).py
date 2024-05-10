import numpy as np
import pyvista as pv
import plotly.graph_objects as go

def get_the_slice(x,y,z, surfacecolor):
    return go.Surface(x=x,
                      y=y,
                      z=z,
                      surfacecolor=surfacecolor,
                      coloraxis='coloraxis')

def get_lims_colors(surfacecolor):# color limits for a slice
    return np.min(surfacecolor), np.max(surfacecolor)


filename = 'plane.vtk'
mesh = pv.read(filename)

mesh.set_active_scalars('PRESSURE_COEFFICIENT')

if hasattr(mesh, 'cells'):
    cells = mesh.cells.reshape(-1, 4)
    is_triangle = cells[:, 0] == 3
    triangles = cells[is_triangle, 1:4]

    pressure_coefficient = mesh.point_data['PRESSURE_COEFFICIENT']
    x=mesh.points[:, 0]
    y=mesh.points[:, 1]
    z=mesh.points[:, 2]

sminz, smaxz = get_lims_colors(pressure_coefficient)
slice_z = get_the_slice(x, y, 0.0 * np.ones(x.shape), pressure_coefficient)

sminy, smaxy = get_lims_colors(pressure_coefficient)
vmin = min([sminz, sminy])
vmax = max([smaxz, smaxy])
slice_y = get_the_slice(x, 0.5 * np.ones(x.shape), z, pressure_coefficient)

def colorax(vmin, vmax):
    return dict(cmin=vmin,
                cmax=vmax)

fig1 = go.Figure(data=[slice_z, slice_y])
fig1.update_layout(
         title_text='Slices in volumetric data', 
         title_x=0.5,
         width=700,
         height=700,
         scene_zaxis_range=[-2,2], 
         coloraxis=dict(colorscale='BrBG',
                        colorbar_thickness=25,
                        colorbar_len=0.75,
                        **colorax(vmin, vmax)))            
      
fig1.show()