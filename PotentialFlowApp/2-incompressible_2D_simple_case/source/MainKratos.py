import KratosMultiphysics
import matplotlib as plt
import numpy as np
from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis import PotentialFlowAnalysis

with open("ProjectParameters.json",'r') as parameter_file:
	parameters = KratosMultiphysics.Parameters(parameter_file.read())

model = KratosMultiphysics.Model()
simulation = PotentialFlowAnalysis(model,parameters)
simulation.Run()


modelpart = model["FluidModelPart.Body2D_Body"]
x = np.zeros(modelpart.NumberOfNodes())
y = np.zeros(modelpart.NumberOfNodes())
z = np.zeros(modelpart.NumberOfNodes())
cp = np.zeros(modelpart.NumberOfNodes())
rho = np.zeros(modelpart.NumberOfNodes())
for i,node in enumerate(modelpart.Nodes):
    x[i] = node.X0 ; y[i] = node.Y0 ; z[i] = node.Z0
    cp[i] = node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
    
# Plot cp vs x
fig,ax  = plt.subplots()
fig.set_figwidth(15.0)
fig.set_figheight(10.0)
ax.plot( x, cp, "o", markersize = 3.0)
ax.grid()
plt.ylabel('Cp')
plt.xlabel('x')
plt.title('Cp vs x')
ax.invert_yaxis()
plt.tight_layout()
fig.savefig("Airfoils_Cp_x.png")
plt.close()