import KratosMultiphysics
from KratosMultiphysics.MeshMovingApplication.mesh_moving_analysis import MeshMovingAnalysis

if __name__ == "__main__":

    with open("ProjectParametersMeshMoving.json", 'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    global_model = KratosMultiphysics.Model()
    mesh_simulation = MeshMovingAnalysis(global_model, parameters)
    mesh_simulation.Run()

main_model_part = global_model.GetModelPart("MainModelPart")
KratosMultiphysics.ModelPartIO("naca0012_3D_0aoa_01", KratosMultiphysics.IO.WRITE | KratosMultiphysics.IO.MESH_ONLY).WriteModelPart(main_model_part)
