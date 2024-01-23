import errno
import numpy as np
import pyvista
import scipy
from matplotlib import pyplot as plt
from scipy.interpolate import RBFInterpolator
import matplotlib.pyplot as plt

def EvaluateTemperature(pos):
    pos_rescaled = RescaleParameter(pos)
    pos_rescaled = pos_rescaled.reshape((1,2))
    q_interpolated = t_interpolator(pos_rescaled).T #note that i need to transpose it back
    t = Phi_t@q_interpolated
    return t

def EvaluateVelocity(pos):
    pos_rescaled = RescaleParameter(pos)
    pos_rescaled = pos_rescaled.reshape((1,2))
    q_interpolated = v_interpolator(pos_rescaled).T #note that i need to transpose it back
    v = Phi_v@q_interpolated
    v = np.reshape(v,(int(v.shape[0]/2),2))
    np.c_[v,np.zeros(v.shape[0])]
    return v

def RescaleParameter(x):
    param_lower_bound = np.array([0.16,3400.0])
    param_upper_bound = np.array([0.2,3900.0])
    for i in range(2):
        y = x.copy() - param_lower_bound
        y = y*1.0/(param_upper_bound-param_lower_bound)
    return y

def TemperatureSlider(value):
    pos[1] = value
    t = EvaluateTemperature(pos)
    mesh.point_data['temperatures'] = t
    v = EvaluateVelocity(pos)
    mesh.point_data['velocities'] = v
    pl.render()

def MassInflowSlider(value):
    pos[0] = value
    t = EvaluateTemperature(pos)
    mesh.point_data['temperatures'] = t
    v = EvaluateVelocity(pos)
    mesh.point_data['velocities'] = v
    pl.render()

with open('simulation_modes.npy', 'rb') as f:
    params = np.load(f)
    Phi_v = np.load(f )
    q_v = np.load(f )
    Phi_t = np.load(f )
    q_t = np.load(f )

with open('simulation_modes.npy', 'rb') as f:
    params = np.load(f)
    Phi_v = np.load(f )
    q_v = np.load(f )
    Phi_t = np.load(f )
    q_t = np.load(f )

plt.scatter(params[0,:],params[1,:])
plt.xlabel('combustion T')
plt.ylabel('mass inflow')
plt.show()

scaled_params = np.zeros(params.shape)
for j in range(params.shape[1]):
    scaled_params[:,j] = RescaleParameter(params[:,j])

t_interpolator = RBFInterpolator(scaled_params.T, q_t.T)
v_interpolator = RBFInterpolator(scaled_params.T, q_v.T)

pos = np.array([0.18, 3350.0])
t = EvaluateTemperature(pos)
v = EvaluateVelocity(pos)

with open('visualization_model.npy', 'rb') as f:
    coords = np.load(f)
    connectivity = np.load(f )

mesh = pyvista.PolyData(coords, connectivity)
mesh.point_data['temperatures'] = t
mesh.point_data['velocities'] = v

pl = pyvista.Plotter()
number = 10
my_cmap = plt.get_cmap('coolwarm',number)
actor = pl.add_mesh(mesh, cmap=my_cmap, scalars='velocities',show_scalar_bar=False)

scalar_bar_actor = pl.add_scalar_bar(title='Original Title', mapper=actor.GetMapper())

pl.add_slider_widget(TemperatureSlider,
    [3400.0,3900.0],
    pointa=(0.025, 0.99),
    pointb=(0.31, 0.99),
    title='Temperature [K]',
    style='modern')

pl.add_slider_widget(MassInflowSlider,
    [0.16,0.20],
    pointa=(0.35, 0.99),
    pointb=(0.64, 0.99),
    title='TotalMassGeneration [kg/s]',
    style='modern')

pl.add_text(
    'Toggle Results',
    position=[70,25], #'lower_left',
    color='blue',
    shadow=False,
    font_size=10,
)

def toggle_results(flag):
    if flag:
        scalar_bar_actor.SetTitle('Velocity [mm/s]')
        max = np.max(mesh.point_data['velocities'])
        pl.update_scalar_bar_range([0,max])
        mesh.set_active_scalars('velocities')
    else:
        scalar_bar_actor.SetTitle('Temperature [K]')
        max = np.max(mesh.point_data['temperatures'])
        pl.update_scalar_bar_range([0,max])
        mesh.set_active_scalars('temperatures')


#toggle_results(False)
_ = pl.add_checkbox_button_widget(toggle_results, value=True)
pl.show(cpos='xy')