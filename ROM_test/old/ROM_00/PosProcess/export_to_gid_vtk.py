import KratosMultiphysics
import KratosMultiphysics.MedApplication as KratosMED
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp 
from KratosMultiphysics.gid_output_process import GiDOutputProcess
from KratosMultiphysics.vtk_output_process import VtkOutputProcess
import numpy as np


def ExportToGiD(model_part,i):
    gid_output = GiDOutputProcess(
        model_part,
            "Gid_errors/output_"+str(i),
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
                "nodal_nonhistorical_results": []
            }
        }"""))
    gid_output.ExecuteInitialize()
    gid_output.ExecuteBeforeSolutionLoop()
    gid_output.ExecuteInitializeSolutionStep()
    gid_output.PrintOutput()
    gid_output.ExecuteFinalizeSolutionStep()
    gid_output.ExecuteFinalize()

def ExportToVtk(model_part,i):
    parameters = KratosMultiphysics.Parameters("""{
                    "model_part_name"                    : "MainModelPart",
                    "output_control_type"                : "step",
                    "output_interval"                    : 1,
                    "file_format"                        : "ascii",
                    "output_precision"                   : 7,
                    "output_sub_model_parts"             : false,
                    "output_path"                        : "vtk_output",
                    "save_output_files_in_folder"        : true,
                    "nodal_solution_step_data_variables" : ["VELOCITY_POTENTIAL","AUXILIARY_VELOCITY_POTENTIAL"],
                    "nodal_data_value_variables"         : []
                    }
                    """)
    parameters["output_path"].SetString("Vtk_errors/vtk_output_"+str(i))

    output = VtkOutputProcess(model_part.GetModel(),parameters)

    output.ExecuteInitialize()
    output.ExecuteBeforeSolutionLoop()
    output.ExecuteInitializeSolutionStep()
    output.PrintOutput()
    output.ExecuteFinalizeSolutionStep()
    output.ExecuteFinalize()


if __name__ == '__main__':

    data_set = np.load("1.8305035058597887, 0.7392086079723365.npy")
    if data_set.ndim != 1:
        if data_set.shape[1] > data_set.shape[0]: data_set = data_set.T
    print(data_set.shape)
    
    if data_set.ndim == 1:
        steps = 1
    else:
        steps = data_set.shape[1]

    for i in range(steps):

        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("MainModelPart")

        model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 2
        model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.0)
        model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 1.0)
        model_part.ProcessInfo.SetValue(KratosMultiphysics.STEP, 1)

        model_part.AddNodalSolutionStepVariable(CPFApp.VELOCITY_POTENTIAL)
        model_part.AddNodalSolutionStepVariable(CPFApp.AUXILIARY_VELOCITY_POTENTIAL)

        # Read the med file
        KratosMED.MedModelPartIO("../salome_mesh_files/model_mesh_0.med", KratosMED.MedModelPartIO.READ).ReadModelPart(model_part)# apply the elements and conditions
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

            if data_set.ndim == 1:
                node.SetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL, data_set[offset])
                node.SetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL, data_set[offset+1])
            else:
                node.SetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL, data_set[offset,i])
                node.SetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL, data_set[offset+1,i])

        ExportToGiD(model_part,i)
        ExportToVtk(model_part,i)