import KratosMultiphysics
import KratosMultiphysics.MedApplication as KratosMED
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp 
from KratosMultiphysics.gid_output_process import GiDOutputProcess
import numpy as np


def ExportToGiD(model_part, mu):
    gid_output = GiDOutputProcess(
        model_part,
            f'RBF_Results/{mu[0]}, {mu[1]}',
        KratosMultiphysics.Parameters("""{
            "result_file_configuration": {
                "gidpost_flags": {
                    "GiDPostMode": "GiD_PostBinary",
                    "WriteDeformedMeshFlag": "WriteUndeformed",
                    "WriteConditionsFlag": "WriteConditions",
                    "MultiFileFlag": "SingleFile"
                },
                "file_label": "step",
                "output_control_type": "step",
                "body_output": true,
                "nodal_results": ["VELOCITY_POTENTIAL","AUXILIARY_VELOCITY_POTENTIAL"],
                "nodal_nonhistorical_results": ["PRESSURE_COEFFICIENT"]
            }
        }"""))
    gid_output.ExecuteInitialize()
    gid_output.ExecuteBeforeSolutionLoop()
    gid_output.ExecuteInitializeSolutionStep()
    gid_output.PrintOutput()
    gid_output.ExecuteFinalizeSolutionStep()
    gid_output.ExecuteFinalize()


if __name__ == '__main__':

    mu = [1.0074408121168927, 0.7093020407269548]

    rbf_data_filename = f"RBF_interpolated_data/{mu[0]}, {mu[1]}.npy"

    data_set = np.load(rbf_data_filename).T

    current_model = KratosMultiphysics.Model()
    model_part = current_model.CreateModelPart("MainModelPart")

    model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 2
    model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.0)
    model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 1.0)
    model_part.ProcessInfo.SetValue(KratosMultiphysics.STEP, 1)

    model_part.AddNodalSolutionStepVariable(CPFApp.VELOCITY_POTENTIAL)
    model_part.AddNodalSolutionStepVariable(CPFApp.AUXILIARY_VELOCITY_POTENTIAL)

    # Read the med file
    KratosMED.MedModelPartIO("model_mesh_0.med", KratosMED.MedModelPartIO.READ).ReadModelPart(model_part)# apply the elements and conditions
    params = KratosMultiphysics.Parameters("""{
        "elements_list"   : [{
            "model_part_name" : "MainModelPart.Parts_Parts_Auto1",
            "element_name"    : "Element2D3N"
        }],
        "conditions_list" : [{
            "model_part_name" : "MainModelPart.PotentialWallCondition2D_Far_field_Auto1",
            "condition_name"  : "WallCondition2D2N"
        },{
            "model_part_name" : "MainModelPart.Body2D_Body",
            "condition_name"  : "WallCondition2D2N"
        }]
    }""")

    modeler = KratosMultiphysics.CreateEntitiesFromGeometriesModeler(current_model, params)
    modeler.SetupModelPart()
    # Assign a fresh properties container to the model
    properties = model_part.CreateNewProperties(1)
    for cond in model_part.Conditions:
        cond.Properties = properties

    for elem in model_part.Elements:
        elem.Properties = properties

    for node in model_part.Nodes:
        offset = np.where(np.arange(1,model_part.NumberOfNodes()+1, dtype=int) == node.Id)[0][0]*2

        node.SetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL, data_set[offset])
        node.SetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL, data_set[offset+1])

    ExportToGiD(model_part, mu)