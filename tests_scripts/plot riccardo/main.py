import numpy as np
import mount_database
import KratosMultiphysics
import pyvista

path = 'C:\\Users\\rross\\Programming\\rocket\\data'
params, velocities, temperatures = mount_database.MountDatabase(path)
# print("params",params.shape,params)
# print("velocities",velocities.shape,velocities)
# print("temperatures",temperatures.shape,temperatures)
print("params.shape=",params.shape)
print("velocities.shape=",velocities.shape)

U,S,V = np.linalg.svd(velocities, full_matrices=False)
print(S)

modes_to_keep = 12
Phi_v = U[:,0:modes_to_keep].copy()
q_v = Phi_v.T@velocities

U,S,V = np.linalg.svd(temperatures, full_matrices=False)
print(S)
Phi_t = U[:,0:modes_to_keep].copy()
q_t = Phi_t.T@temperatures

#########################
with open('simulation_modes.npy', 'wb') as f:
    np.save(f, params)
    np.save(f, Phi_v)
    np.save(f, q_v)
    np.save(f, Phi_t)
    np.save(f, q_t)

# print("reduced repr of vel", Phi_v.T@velocities)



#plotting velocity on step i
i=0
print(velocities)

v = velocities[:,i].copy()
v = np.reshape(v,(int(v.shape[0]/2),2))
np.c_[v,np.zeros(v.shape[0])]

t = temperatures[:,i]
print(t)


print(dir(KratosMultiphysics.VariableUtils()))

#######################################################
current_model = KratosMultiphysics.Model()
all_mp= current_model.CreateModelPart("Main")
all_mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
all_mp.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
all_mp.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)

model_part_io = KratosMultiphysics.ModelPartIO("rocket_ext_edit")
model_part_io.ReadModelPart(all_mp)

mp = all_mp.GetSubModelPart("Rocket")

##nasty, but i need to do a reorder, so that nodes in the submodelpart are numbered first
for node in all_mp.Nodes:
    node.SetValue(KratosMultiphysics.AUX_INDEX,-1)

counter = 0
for node in mp.Nodes:
    node.SetValue(KratosMultiphysics.AUX_INDEX,counter)
    counter+=1
for node in all_mp.Nodes:
    aux_id = node.GetValue(KratosMultiphysics.AUX_INDEX)
    if(aux_id < 0):
        node.SetValue(KratosMultiphysics.AUX_INDEX,counter)
        counter+=1

print(mp)
coords = np.array(KratosMultiphysics.VariableUtils().GetCurrentPositionsVector(mp.Nodes,3))
nnodes = len(mp.Nodes)
coords = np.reshape(coords,(nnodes,3))

connectivity = []
for elem in mp.Elements:
    geom = elem.GetGeometry()
    nn = len(geom)
    connectivity.append(nn)
    for node in geom:
        connectivity.append(int(node.GetValue(KratosMultiphysics.AUX_INDEX)))
connectivity=np.array(connectivity)

with open('visualization_model.npy', 'wb') as f:
    np.save(f, coords)
    np.save(f, connectivity)

print("aaaaaaaaaa")

mesh = pyvista.PolyData(coords, connectivity)
# mesh.point_data.set_array(v,'velocities')
# mesh.point_data.set_array(t,'temperatures')
mesh.point_data['temperatures'] = t
mesh.point_data['velocities'] = v

mesh.set_active_scalars('temperatures')

pl = pyvista.Plotter()
pl.add_mesh(mesh, cmap='coolwarm')

def TemperatureSlider(value):
    res = value
    return

def MassInflowSlider(value):
    res = value
    return

pl.add_slider_widget(TemperatureSlider,
    [3300.0,3500.0],
    pointa=(0.025, 0.99),
    pointb=(0.31, 0.99),
    title='Temperature [K]',
    style='modern')
pl.add_slider_widget(MassInflowSlider,
    [0.16,0.18],
    pointa=(0.35, 0.99),
    pointb=(0.64, 0.99),
    title='TotalMassGeneration [kg/s]',
    style='modern')

def show_temperatures(flag):
    print("1111111111111")
    print(pl)
    mesh.set_active_scalars('temperatures')
    # pl.update()
    pl.update_scalars("temperatures", mesh=mesh)
    print(pl)
    print("22222222222222")

def show_velocities(flag):

    mesh.set_active_scalars('velocities')


_ = pl.add_checkbox_button_widget(show_velocities, value=True)
print("11111")
pl.show(cpos='xy', interactive=True) #, show_edges=False)
#interactive=True,
#, interactive_update=True
print("bbbbbbb")

