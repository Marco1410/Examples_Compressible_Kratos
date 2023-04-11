import KratosMultiphysics
from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis import PotentialFlowAnalysis

with open("ProjectParametersPrimal.json",'r') as parameter_file:
	parameters = KratosMultiphysics.Parameters(parameter_file.read())

model = KratosMultiphysics.Model()
simulation = PotentialFlowAnalysis(model,parameters)
simulation.Run()
