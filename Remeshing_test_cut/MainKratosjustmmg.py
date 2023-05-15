import KratosMultiphysics
import KratosMultiphysics.MeshingApplication

def PrintOutput(main_model_part,filename):
    from KratosMultiphysics.gid_output_process import GiDOutputProcess
    gid_output = GiDOutputProcess(main_model_part,
                                filename,
                                KratosMultiphysics.Parameters("""
                                    {
                                        "result_file_configuration" : {
                                            "gidpost_flags": {
                                                "GiDPostMode": "GiD_PostBinary",
                                                "WriteDeformedMeshFlag": "WriteUndeformed",
                                                "WriteConditionsFlag": "WriteElementsOnly",
                                                "MultiFileFlag": "SingleFile"
                                            },
                                            "nodal_results": ["DISTANCE"],
                                            "nodal_nonhistorical_results": ["NODAL_H"]
                                        }
                                    }
                                    """)
                                )

    gid_output.ExecuteInitialize()
    gid_output.ExecuteBeforeSolutionLoop()
    gid_output.ExecuteInitializeSolutionStep()
    gid_output.PrintOutput()
    gid_output.ExecuteFinalizeSolutionStep()
    gid_output.ExecuteFinalize()

def PrintOutputMetric(main_model_part,filename):

    for node in main_model_part.Nodes:
        tensor_3d = node.GetValue(KratosMultiphysics.MeshingApplication.METRIC_TENSOR_3D)
        node.SetValue(KratosMultiphysics.TEMPERATURE, tensor_3d.norm_2())
    from KratosMultiphysics.gid_output_process import GiDOutputProcess
    gid_output = GiDOutputProcess(main_model_part,
                                filename,
                                KratosMultiphysics.Parameters("""
                                    {
                                        "result_file_configuration" : {
                                            "gidpost_flags": {
                                                "GiDPostMode": "GiD_PostBinary",
                                                "WriteDeformedMeshFlag": "WriteUndeformed",
                                                "WriteConditionsFlag": "WriteElementsOnly",
                                                "MultiFileFlag": "SingleFile"
                                            },
                                            "nodal_results": ["DISTANCE"],
                                            "nodal_nonhistorical_results": ["NODAL_H", "TEMPERATURE","DISTANCE_GRADIENT"]
                                        }
                                    }
                                    """)
                                )

    gid_output.ExecuteInitialize()
    gid_output.ExecuteBeforeSolutionLoop()
    gid_output.ExecuteInitializeSolutionStep()
    gid_output.PrintOutput()
    gid_output.ExecuteFinalizeSolutionStep()
    gid_output.ExecuteFinalize()


def _CalculateDistance(main_model_part, skin_model_part):

        calculate_distance_process = KratosMultiphysics.CalculateDistanceToSkinProcess3D(
                main_model_part,
                skin_model_part)
        calculate_distance_process.Execute()

def _ComputeLevelSetMetric(main_model_part,min_size_level):
        # Extending distace to all the domain
        import KratosMultiphysics.python_linear_solver_factory #Linear solver for variational distance process
        linear_solver_settings=KratosMultiphysics.Parameters("""
        {
            "solver_type": "amgcl",
            "max_iteration": 400,
            "gmres_krylov_space_dimension": 1000,
            "smoother_type":"ilu0",
            "coarsening_type":"ruge_stuben",
            "coarse_enough" : 5000,
            "krylov_type": "lgmres",
            "tolerance": 1e-9,
            "verbosity": 0,
            "scaling": false
        }""")
        linear_solver = KratosMultiphysics.python_linear_solver_factory.ConstructSolver(linear_solver_settings)
        maximum_iterations = 2
        variational_distance_process = KratosMultiphysics.VariationalDistanceCalculationProcess3D(
            main_model_part,
            linear_solver,
            maximum_iterations)
        variational_distance_process.Execute()

        #COMPUTE DISTANCE GRADIENT AND NODAL_H
        local_gradient = KratosMultiphysics.ComputeNodalGradientProcess3D(main_model_part,
            KratosMultiphysics.DISTANCE,
            KratosMultiphysics.DISTANCE_GRADIENT,
            KratosMultiphysics.NODAL_AREA)
        local_gradient.Execute()

        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(main_model_part)
        find_nodal_h.Execute()
        PrintOutput(main_model_part,'extended_distance_mesh')
                # size distribution, matrix of sizes [level_set, size]
        metric_parameters = KratosMultiphysics.Parameters("""
        {
            "sizing_parameters": {

                "reference_variable_name"               : "DISTANCE",
                "boundary_layer_max_distance"           : 60,
                "interpolation"                         : "piecewise_linear",
                "size_distribution": [
                                                [ -60,100  ],
                                                [-0.5,0.1],
                                                [   0,0.005],
                                                [ 0.5,0.1],
                                                [  60,100  ]
                                     ]
            },
            "enforce_current"                      : false,
            "anisotropy_remeshing"                 : false

        }
        """)

        # Calculate level set metric to refine according to geometry level-set
        metric_process = KratosMultiphysics.MeshingApplication.ComputeLevelSetSolMetricProcess3D(main_model_part,  KratosMultiphysics.DISTANCE_GRADIENT, metric_parameters)
        metric_process.Execute()

def _RemeshAfterCut(main_model_part, skin_model_part):
        # Perform refinement before isosurface cutting
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(KratosMultiphysics.MeshingApplication.METRIC_TENSOR_3D,main_model_part.Nodes)
        _CalculateDistance(main_model_part, skin_model_part)
        _ComputeLevelSetMetric(main_model_part)
        PrintOutput(main_model_part,'extended_distance_after_cut_mesh')

        mmg_parameters = KratosMultiphysics.Parameters("""
        {
            "discretization_type"              : "STANDARD",
            "save_external_files"              : false,
            "interpolate_nodal_values"         : false,
            "preserve_flags"                   : false,
            "initialize_entities"              : false,
            "echo_level"                       : 5
        }
        """)
        mmg_process = KratosMultiphysics.MeshingApplication.MmgProcess3D(main_model_part, mmg_parameters)
        mmg_process.Execute()
        PrintOutput(main_model_part,'remeshed_levelset_output_after_cut')



###
### 1. skin and volume mdpa is read.
###
model=KratosMultiphysics.Model()

# # Read skin file  (only triangles)
file_name = "naca0012_3D_skin.gid/naca0012_3D_skin"
skin_model_part=model.CreateModelPart("Skin")
KratosMultiphysics.ModelPartIO(file_name).ReadModelPart(skin_model_part)

#file_name="../NavierStokes/naca0012_3D.gid/naca0012_3D"
file_name= "naca0012_3D_farfield.gid/naca0012_3D_farfield"
model=KratosMultiphysics.Model()
main_model_part=model.CreateModelPart("MainModelPart")
main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3
main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)
main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FLAG_VARIABLE)
KratosMultiphysics.ModelPartIO(file_name,KratosMultiphysics.IO.READ| KratosMultiphysics.IO.MESH_ONLY).ReadModelPart(main_model_part)

min_size_level=0.5

###
### 2. We compute  the metric according to the level set (skin).
###
KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(KratosMultiphysics.MeshingApplication.METRIC_TENSOR_3D,main_model_part.Nodes)
_CalculateDistance(main_model_part, skin_model_part)
_ComputeLevelSetMetric(main_model_part,min_size_level)


###
### 3. Perform isosurface cutting
###

mmg_parameters = KratosMultiphysics.Parameters("""
{
    "echo_level"                       : 3,
    "discretization_type"              : "IsoSurface",
    "interpolate_nodal_values"         : false,
    "save_external_files"              : false,
    "initialize_entities"              : false,
    "preserve_flags"                   : false,
    "isosurface_parameters"            :
    {
            "isosurface_variable"              : "DISTANCE",
            "nonhistorical_variable"           : false,
            "use_metric_field"                 : true,
            "remove_internal_regions"          : true
    },
    "force_sizes"                      :
    {
        "force_min"                           : true,
        "minimal_size"                        : 0.005,
        "force_max"                           : true,
        "maximal_size"                        : 100
    }


}
""")

mmg_process = KratosMultiphysics.MeshingApplication.MmgProcess3D(main_model_part, mmg_parameters)
mmg_process.Execute()
PrintOutputMetric(main_model_part,"output")
tmoc = KratosMultiphysics.TetrahedralMeshOrientationCheck
throw_errors = False
flags = (tmoc.COMPUTE_NODAL_NORMALS).AsFalse() | (tmoc.COMPUTE_CONDITION_NORMALS).AsFalse()
flags |= tmoc.ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS
KratosMultiphysics.TetrahedralMeshOrientationCheck(main_model_part,throw_errors, flags).Execute()

PrintOutput(main_model_part,'remeshed_output')
KratosMultiphysics.ModelPartIO("remeshed_final_cut_mesh", KratosMultiphysics.IO.WRITE | KratosMultiphysics.IO.MESH_ONLY).WriteModelPart(main_model_part)



# TODO: It works for an old version of mmg, now there are some problems with it, so be carefull with the corners
#  "advanced_parameters"                  :
#     {

#         "local_entity_parameters_list"        : [
#             {
#                 "model_part_name_list" : ["Wall2", "Inlet", "Outlet", "Wall1","Top"],
#                 "hmin"            : 15,
#                 "hmax"            : 20,
#                 "hausdorff_value" : 10000000.0
#             }],
#         "deactivate_detect_angle"             : false,
#         "force_gradation_value"               : false,
#         "mesh_optimization_only"              : false,
#         "gradation_value"                     : 1.9
#     }
