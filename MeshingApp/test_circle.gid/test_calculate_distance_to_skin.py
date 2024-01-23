import KratosMultiphysics
import KratosMultiphysics.MeshingApplication
import KratosMultiphysics.python_linear_solver_factory
from KratosMultiphysics.gid_output_process import GiDOutputProcess

def test_circle_calculate_distance_to_skin_2d():
        
    # Set the problem domain using the structured mesh generator process
    current_model = KratosMultiphysics.Model()
    model_part = current_model.CreateModelPart("ModelPart")

    model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 2
    model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.0)
    model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 1.0)
    model_part.ProcessInfo.SetValue(KratosMultiphysics.STEP, 1)

    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)

    problem_domain = KratosMultiphysics.Quadrilateral2D4(
        KratosMultiphysics.Node(1, 0.0, 0.0, 0.0),
        KratosMultiphysics.Node(2, 0.0, 1.0, 0.0),
        KratosMultiphysics.Node(3, 1.0, 1.0, 0.0),
        KratosMultiphysics.Node(4, 1.0, 0.0, 0.0))
    parameters = KratosMultiphysics.Parameters("{}")
    parameters.AddEmptyValue("element_name").SetString("Element2D3N")
    parameters.AddEmptyValue("condition_name").SetString("LineCondition2D2N")
    parameters.AddEmptyValue("create_skin_sub_model_part").SetBool(False)
    parameters.AddEmptyValue("number_of_divisions").SetInt(50)
    KratosMultiphysics.StructuredMeshGeneratorProcess(problem_domain, model_part, parameters).Execute()
        
    # Set geometry
    skin_model_part = current_model.CreateModelPart("circle")
    KratosMultiphysics.ModelPartIO("test_circle.gid/test_circle").ReadModelPart(skin_model_part)
    
    # Call the CalculateDistanceToSkinProcess()
    KratosMultiphysics.CalculateDistanceToSkinProcess2D(model_part, skin_model_part).Execute()
        
    # Extending distace to all the domain
    linear_solver_settings=KratosMultiphysics.Parameters("""
    {   
        "solver_type": "LinearSolversApplication.pardiso_lu"
    }""")
    linear_solver = KratosMultiphysics.python_linear_solver_factory.ConstructSolver(linear_solver_settings)
    maximum_iterations = 1
    utility = KratosMultiphysics.VariationalDistanceCalculationProcess2D(
        model_part,
        linear_solver,
        maximum_iterations)
    utility.Execute()
        
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
    KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(KratosMultiphysics.MeshingApplication.METRIC_TENSOR_2D, model_part.Nodes)
    # size distribution, matrix of sizes [level_set, size]
    metric_parameters = KratosMultiphysics.Parameters("""
    {
        "minimal_size"                         : 0.001,
        "maximal_size"                         : 1.0,
            "sizing_parameters": {
                "reference_variable_name"               : "DISTANCE",
                "boundary_layer_max_distance"           : 1,
                "interpolation"                         : "linear"
                },
                "enforce_current"                      : false,
                "anisotropy_remeshing"                 : false
    }
    """)
    # Calculate level set metric to refine according to geometry level-set
    metric_process = KratosMultiphysics.MeshingApplication.ComputeLevelSetSolMetricProcess2D(model_part,  KratosMultiphysics.DISTANCE_GRADIENT, metric_parameters)
    metric_process.Execute()

    for node in model_part.Nodes:
        tensor_3d = node.GetValue(KratosMultiphysics.MeshingApplication.METRIC_TENSOR_2D)
        node.SetValue(KratosMultiphysics.TEMPERATURE, tensor_3d.norm_2())
            
    gid_output = GiDOutputProcess(
        model_part,
            "output_metric",
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

    mmg_parameters = KratosMultiphysics.Parameters("""
        {
            "discretization_type": "Isosurface",
            "force_sizes": {
                "force_max": true,
                "force_min": true,
                "maximal_size": 1.0,
                "minimal_size": 0.001
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
    mmg_process = KratosMultiphysics.MeshingApplication.MmgProcess2D(model_part, mmg_parameters)
    # We remesh
    mmg_process.Execute()

    gid_output = GiDOutputProcess(
        model_part,
            "final_output",
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

def test_sphere_calculate_distance_to_skin_3d():

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
        KratosMultiphysics.Node(1, -0.50, -0.75, -0.22),
        KratosMultiphysics.Node(2,  1.50, -0.75, -0.22),
        KratosMultiphysics.Node(3,  1.50,  0.75, -0.22),
        KratosMultiphysics.Node(4, -0.50,  0.75, -0.22),
        KratosMultiphysics.Node(5, -0.50, -0.75,  0.22),
        KratosMultiphysics.Node(6,  1.50, -0.75,  0.22),
        KratosMultiphysics.Node(7,  1.50,  0.75,  0.22),
        KratosMultiphysics.Node(8, -0.50,  0.75,  0.22))
    parameters = KratosMultiphysics.Parameters("{}")
    parameters.AddEmptyValue("element_name").SetString("Element3D4N")
    parameters.AddEmptyValue("condition_name").SetString("SurfaceCondition3D3N")
    parameters.AddEmptyValue("create_skin_sub_model_part").SetBool(False)
    parameters.AddEmptyValue("number_of_divisions").SetInt(50)
    KratosMultiphysics.StructuredMeshGeneratorProcess(problem_domain, model_part, parameters).Execute()
    
    # Set aerofoil geometry
    skin_model_part = current_model.CreateModelPart("Sphere")
    KratosMultiphysics.ModelPartIO("test_sphere").ReadModelPart(skin_model_part)

    iterations = 1
    i = 0

    for i in range(iterations):
        i += 1
    
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
            "maximal_size"                         : 1.0,
            "sizing_parameters": {
                "reference_variable_name"               : "DISTANCE",
                "boundary_layer_max_distance"           : 1,
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
                "maximal_size": 1.0,
                "minimal_size": 0.005
            }
        }
        """)

        # We create the remeshing utility
        mmg_process = KratosMultiphysics.MeshingApplication.MmgProcess3D(model_part, mmg_parameters)
        # We remesh
        mmg_process.Execute()

        tmoc = KratosMultiphysics.TetrahedralMeshOrientationCheck
        throw_errors = True
        flags = (tmoc.COMPUTE_NODAL_NORMALS).AsFalse() | (tmoc.COMPUTE_CONDITION_NORMALS).AsFalse()
        flags |= tmoc.ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS
        KratosMultiphysics.TetrahedralMeshOrientationCheck(model_part,throw_errors, flags).Execute()

        for node in model_part.Nodes:
            tensor_3d = node.GetValue(KratosMultiphysics.MeshingApplication.METRIC_TENSOR_3D)
            node.SetValue(KratosMultiphysics.TEMPERATURE, tensor_3d.norm_2())

        gid_output = GiDOutputProcess(
            model_part,
            "output"+str(i),
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
        "maximal_size"                         : 1.0,
        "sizing_parameters": {
            "reference_variable_name"               : "DISTANCE",
            "boundary_layer_max_distance"           : 1,
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
                "maximal_size": 1.0,
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

    tmoc = KratosMultiphysics.TetrahedralMeshOrientationCheck
    throw_errors = True
    flags = (tmoc.COMPUTE_NODAL_NORMALS).AsFalse() | (tmoc.COMPUTE_CONDITION_NORMALS).AsFalse()
    flags |= tmoc.ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS
    KratosMultiphysics.TetrahedralMeshOrientationCheck(model_part,throw_errors, flags).Execute()

    for node in model_part.Nodes:
        tensor_3d = node.GetValue(KratosMultiphysics.MeshingApplication.METRIC_TENSOR_3D)
        node.SetValue(KratosMultiphysics.TEMPERATURE, tensor_3d.norm_2())

    # Print results (left it here for debugging)
    gid_output = GiDOutputProcess(
        model_part,
            "final_output",
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

if __name__ == '__main__':
    domain_size = 2
    if (domain_size == 2):
        test_circle_calculate_distance_to_skin_2d()
    else:
        test_sphere_calculate_distance_to_skin_3d()
