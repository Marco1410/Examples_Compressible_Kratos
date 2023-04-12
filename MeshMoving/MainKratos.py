import KratosMultiphysics
from KratosMultiphysics.MeshMovingApplication.mesh_moving_analysis import MeshMovingAnalysis

with open("ProjectParametersPrimal.json",'r') as parameter_file:
	parameters = KratosMultiphysics.Parameters(parameter_file.read())

model = KratosMultiphysics.Model()
simulation = MeshMovingAnalysis(model,parameters)
simulation.Run()