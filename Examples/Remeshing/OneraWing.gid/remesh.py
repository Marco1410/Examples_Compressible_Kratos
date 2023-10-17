import KratosMultiphysics
import KratosMultiphysics.MeshingApplication
import KratosMultiphysics.python_linear_solver_factory
from KratosMultiphysics.gid_output_process import GiDOutputProcess

# Set the problem domain using the structured mesh generator process
current_model = KratosMultiphysics.Model()
model_part = current_model.CreateModelPart("ModelPart")

model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3
model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.0)
model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 1.0)
model_part.ProcessInfo.SetValue(KratosMultiphysics.STEP, 1)

model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)
model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)

problem_domain = KratosMultiphysics.Hexahedra3D8(
    KratosMultiphysics.Node(1, -4.0, 0.0, -3.0),
    KratosMultiphysics.Node(2,  6.0, 0.0, -3.0),
    KratosMultiphysics.Node(3,  6.0, 3.0, -3.0),
    KratosMultiphysics.Node(4, -4.0, 3.0, -3.0),
    KratosMultiphysics.Node(5, -4.0, 0.0,  3.0),
    KratosMultiphysics.Node(6,  6.0, 0.0,  3.0),
    KratosMultiphysics.Node(7,  6.0, 3.0,  3.0),
    KratosMultiphysics.Node(8, -4.0, 3.0,  3.0))
parameters = KratosMultiphysics.Parameters("{}")
parameters.AddEmptyValue("element_name").SetString("Element3D4N")
parameters.AddEmptyValue("condition_name").SetString("SurfaceCondition3D3N")
parameters.AddEmptyValue("create_skin_sub_model_part").SetBool(False)
parameters.AddEmptyValue("number_of_divisions").SetInt(50)
KratosMultiphysics.StructuredMeshGeneratorProcess(problem_domain, model_part, parameters).Execute()

# KratosMultiphysics.ModelPartIO("farfield", KratosMultiphysics.IO.WRITE | KratosMultiphysics.IO.MESH_ONLY).WriteModelPart(model_part)
# KratosMultiphysics.ModelPartIO("farfield").ReadModelPart(model_part)

# Set aerofoil geometry
skin_model_part = current_model.CreateModelPart("Aerofoil")
KratosMultiphysics.ModelPartIO("OneraWing").ReadModelPart(skin_model_part)

iterations = 5
i = 0

for i in range(iterations):
    i += 1
    
    # Call the CalculateDistanceToSkinProcess()
    KratosMultiphysics.CalculateDistanceToSkinProcess3D(model_part, skin_model_part, 1e-12).Execute()

    gid_output = GiDOutputProcess(
        model_part,
        "Results/output_distance"+str(i),
        KratosMultiphysics.Parameters("""{
            "result_file_configuration": {
                "gidpost_flags": {
                    "GiDPostMode": "GiD_PostAscii",
                    "WriteDeformedMeshFlag": "WriteUndeformed",
                    "WriteConditionsFlag": "WriteConditions",
                    "MultiFileFlag": "SingleFile"
                },
                "file_label": "time",
                "output_control_type": "step",
                "body_output": true,
                "nodal_results": ["DISTANCE","DISTANCE_GRADIENT"],
                "nodal_nonhistorical_results": ["NODAL_H", "TEMPERATURE","NODAL_AREA"]
            }
        }"""))
    gid_output.ExecuteInitialize()
    gid_output.ExecuteBeforeSolutionLoop()
    gid_output.ExecuteInitializeSolutionStep()
    gid_output.PrintOutput()
    gid_output.ExecuteFinalizeSolutionStep()
    gid_output.ExecuteFinalize()

    # Extending distace to all the domain
    linear_solver_settings=KratosMultiphysics.Parameters("""
    {   
        "solver_type": "LinearSolversApplication.pardiso_lu"
    }""")
    linear_solver = KratosMultiphysics.python_linear_solver_factory.ConstructSolver(linear_solver_settings)
    maximum_iterations = 1
    KratosMultiphysics.VariationalDistanceCalculationProcess3D(
        model_part,
        linear_solver,
        maximum_iterations).Execute()
 
    gid_output = GiDOutputProcess(
        model_part,
        "Results/output_distance_extended"+str(i),
        KratosMultiphysics.Parameters("""{
            "result_file_configuration": {
                "gidpost_flags": {
                    "GiDPostMode": "GiD_PostAscii",
                    "WriteDeformedMeshFlag": "WriteUndeformed",
                    "WriteConditionsFlag": "WriteConditions",
                    "MultiFileFlag": "SingleFile"
                },
                "file_label": "time",
                "output_control_type": "step",
                "body_output": true,
                "nodal_results": ["DISTANCE","DISTANCE_GRADIENT"],
                "nodal_nonhistorical_results": ["NODAL_H", "TEMPERATURE","NODAL_AREA"]
            }
        }"""))
    gid_output.ExecuteInitialize()
    gid_output.ExecuteBeforeSolutionLoop()
    gid_output.ExecuteInitializeSolutionStep()
    gid_output.PrintOutput()
    gid_output.ExecuteFinalizeSolutionStep()
    gid_output.ExecuteFinalize()

    find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(model_part)
    find_nodal_h.Execute()

    KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(KratosMultiphysics.NODAL_AREA, model_part.Nodes)
    #COMPUTE DISTANCE GRADIENT
    local_gradient = KratosMultiphysics.ComputeNodalGradientProcess3D(model_part,
        KratosMultiphysics.DISTANCE,
        KratosMultiphysics.DISTANCE_GRADIENT,
        KratosMultiphysics.NODAL_AREA)
    local_gradient.Execute()

    # We set to zero the metric
    KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(KratosMultiphysics.MeshingApplication.METRIC_TENSOR_3D, model_part.Nodes)
    # size distribution, matrix of sizes [level_set, size]
    metric_parameters = KratosMultiphysics.Parameters("""
    {
        "minimal_size"                         : 0.005,
        "maximal_size"                         : 10.0,
        "sizing_parameters": {
            "reference_variable_name"               : "DISTANCE",
            "boundary_layer_max_distance"           : 10,
            "interpolation"                         : "linear",
            "size_distribution": [
                                    []
                                 ]
            },
            "enforce_current"                      : false,
            "anisotropy_remeshing"                 : false
    }
    """)
    # Calculate level set metric to refine according to geometry level-set
    metric_process = KratosMultiphysics.MeshingApplication.ComputeLevelSetSolMetricProcess3D(model_part,  KratosMultiphysics.DISTANCE_GRADIENT, metric_parameters)
    metric_process.Execute()

    mmg_parameters = KratosMultiphysics.Parameters("""
    {
        "force_sizes": {
            "force_max": true,
            "force_min": true,
            "maximal_size": 10.0,
            "minimal_size": 0.005
        }
    }
    """)

    # We create the remeshing utility
    mmg_process = KratosMultiphysics.MeshingApplication.MmgProcess3D(model_part, mmg_parameters)
    # We remesh
    mmg_process.Execute()

    for node in model_part.Nodes:
        tensor_3d = node.GetValue(KratosMultiphysics.MeshingApplication.METRIC_TENSOR_3D)
        node.SetValue(KratosMultiphysics.TEMPERATURE, tensor_3d.norm_2())

    gid_output = GiDOutputProcess(
        model_part,
        "Results/output"+str(i),
        KratosMultiphysics.Parameters("""{
            "result_file_configuration": {
                "gidpost_flags": {
                    "GiDPostMode": "GiD_PostAscii",
                    "WriteDeformedMeshFlag": "WriteUndeformed",
                    "WriteConditionsFlag": "WriteConditions",
                    "MultiFileFlag": "SingleFile"
                },
                "file_label": "time",
                "output_control_type": "step",
                "body_output": true,
                "nodal_results": ["DISTANCE","DISTANCE_GRADIENT"],
                "nodal_nonhistorical_results": ["NODAL_H", "TEMPERATURE","NODAL_AREA"]
            }
        }"""))
    gid_output.ExecuteInitialize()
    gid_output.ExecuteBeforeSolutionLoop()
    gid_output.ExecuteInitializeSolutionStep()
    gid_output.PrintOutput()
    gid_output.ExecuteFinalizeSolutionStep()
    gid_output.ExecuteFinalize()

# Call the CalculateDistanceToSkinProcess()
KratosMultiphysics.CalculateDistanceToSkinProcess3D(model_part, skin_model_part, 1e-12).Execute()
        
# Extending distace to all the domain
linear_solver_settings=KratosMultiphysics.Parameters("""
{   
    "solver_type": "LinearSolversApplication.pardiso_lu"
}""")
linear_solver = KratosMultiphysics.python_linear_solver_factory.ConstructSolver(linear_solver_settings)
maximum_iterations = 1
KratosMultiphysics.VariationalDistanceCalculationProcess3D(
    model_part,
    linear_solver,
    maximum_iterations).Execute()
    
find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(model_part)
find_nodal_h.Execute()

KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(KratosMultiphysics.NODAL_AREA, model_part.Nodes)
#COMPUTE DISTANCE GRADIENT
local_gradient = KratosMultiphysics.ComputeNodalGradientProcess2D(model_part,
    KratosMultiphysics.DISTANCE,
    KratosMultiphysics.DISTANCE_GRADIENT,
    KratosMultiphysics.NODAL_AREA)
local_gradient.Execute()

# We set to zero the metric
KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(KratosMultiphysics.MeshingApplication.METRIC_TENSOR_3D, model_part.Nodes)
# size distribution, matrix of sizes [level_set, size]
metric_parameters = KratosMultiphysics.Parameters("""
{
    "minimal_size"                         : 0.005,
    "maximal_size"                         : 10.0,
    "sizing_parameters": {
        "reference_variable_name"               : "DISTANCE",
        "boundary_layer_max_distance"           : 10,
        "interpolation"                         : "linear"
        },
        "enforce_current"                      : false,
        "anisotropy_remeshing"                 : false
}
""")
# Calculate level set metric to refine according to geometry level-set
metric_process = KratosMultiphysics.MeshingApplication.ComputeLevelSetSolMetricProcess3D(model_part,  KratosMultiphysics.DISTANCE_GRADIENT, metric_parameters)
metric_process.Execute()

mmg_parameters = KratosMultiphysics.Parameters("""
    {
        "discretization_type": "Isosurface",
        "force_sizes": {
            "force_max": true,
            "force_min": true,
            "maximal_size": 10.0,
            "minimal_size": 0.005
        },
        "advanced_parameters"                  : {
            "normal_regularization_mesh"       : true,
            "mesh_optimization_only"           : true
        },
        "isosurface_parameters": {
            "invert_value": false,
            "isosurface_variable": "DISTANCE",
            "nonhistorical_variable": false,
            "remove_internal_regions": true,
            "use_metric_field": true
        }
    }
    """)

# We create the remeshing utility
mmg_process = KratosMultiphysics.MeshingApplication.MmgProcess3D(model_part, mmg_parameters)
# We remesh
mmg_process.Execute()

for node in model_part.Nodes:
    tensor_3d = node.GetValue(KratosMultiphysics.MeshingApplication.METRIC_TENSOR_3D)
    node.SetValue(KratosMultiphysics.TEMPERATURE, tensor_3d.norm_2())

gid_output = GiDOutputProcess(
    model_part,
        "Results/final_output",
    KratosMultiphysics.Parameters("""{
        "result_file_configuration": {
            "gidpost_flags": {
                "GiDPostMode": "GiD_PostAscii",
                "WriteDeformedMeshFlag": "WriteUndeformed",
                "WriteConditionsFlag": "WriteConditions",
                "MultiFileFlag": "SingleFile"
            },
            "file_label": "time",
            "output_control_type": "step",
            "body_output": true,
            "nodal_results": ["DISTANCE","DISTANCE_GRADIENT"],
            "nodal_nonhistorical_results": ["NODAL_H", "TEMPERATURE","NODAL_AREA"]
        }
    }"""))
gid_output.ExecuteInitialize()
gid_output.ExecuteBeforeSolutionLoop()
gid_output.ExecuteInitializeSolutionStep()
gid_output.PrintOutput()
gid_output.ExecuteFinalizeSolutionStep()
gid_output.ExecuteFinalize()