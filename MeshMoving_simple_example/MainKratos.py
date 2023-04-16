import KratosMultiphysics
from KratosMultiphysics.MeshMovingApplication.mesh_moving_analysis import MeshMovingAnalysis

with open("ProjectParametersPrimal.json",'r') as parameter_file:
	parameters = KratosMultiphysics.Parameters(parameter_file.read())

model = KratosMultiphysics.Model()
simulation = MeshMovingAnalysis(model,parameters)
simulation.Run()

model_part = model.GetModelPart("MainModelPart")
KratosMultiphysics.ModelPartIO("naca0012_newaoa", KratosMultiphysics.IO.WRITE | KratosMultiphysics.IO.MESH_ONLY).WriteModelPart(model_part)