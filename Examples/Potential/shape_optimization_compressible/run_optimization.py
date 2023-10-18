# Import Kratos core and apps
import KratosMultiphysics
import KratosMultiphysics.ShapeOptimizationApplication.optimizer_factory as optimizer_factory

#Disable logs
KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

with open("optimization_parameters.json",'r') as parameter_file:
    parameters_optimization = KratosMultiphysics.Parameters(parameter_file.read())


model = KratosMultiphysics.Model()
optimizer = optimizer_factory.CreateOptimizer(parameters_optimization["optimization_settings"], model)
optimizer.Optimize()

body_model_part = model.GetModelPart("MainModelPart.Body2D_Body")
KratosMultiphysics.ModelPartIO("final_bodyfitted_model_part", KratosMultiphysics.IO.WRITE | KratosMultiphysics.IO.MESH_ONLY).WriteModelPart(body_model_part)
