import KratosMultiphysics
from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis import PotentialFlowAnalysis

with open("ProjectParametersPrimal.json",'r') as parameter_file:
	parameters = KratosMultiphysics.Parameters(parameter_file.read())

model = KratosMultiphysics.Model()
simulation = PotentialFlowAnalysis(model,parameters)
simulation.Run()

# model_part = model.GetModelPart("MainModelPart")
# KratosMultiphysics.ModelPartIO("final_model_part", KratosMultiphysics.IO.WRITE | KratosMultiphysics.IO.MESH_ONLY).WriteModelPart(model_part)