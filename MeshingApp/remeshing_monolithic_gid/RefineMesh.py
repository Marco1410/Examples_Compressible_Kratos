import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.MeshingApplication as KratosMeshing



def CreateGidControlOutput(output_name, model_part):
    from KratosMultiphysics.gid_output_process import GiDOutputProcess
    gid_output = GiDOutputProcess(
            model_part,
            output_name,
            KratosMultiphysics.Parameters("""
                {
                    "result_file_configuration" : {
                        "gidpost_flags": {
                            "GiDPostMode": "GiD_PostBinary",
                            "MultiFileFlag": "SingleFile"
                        },
                        "nodal_results"       : [],
                        "nodal_nonhistorical_results": ["NODAL_H","METRIC_TENSOR_2D","VELOCITY"],
                        "nodal_flags_results": []
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



model = KratosMultiphysics.Model()
main_model_part = model.CreateModelPart("main")
KratosMultiphysics.ModelPartIO("RectangularCylinder_75k").ReadModelPart(main_model_part)
main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)


with open("average_velocity_field_RectangularCylinder2D_75k_600.0.dat") as dat_file:
    lines=dat_file.readlines()
    for line, node in zip(lines, main_model_part.Nodes):
        velocity = KratosMultiphysics.Vector(3, 0.0)
        velocity[0] = float(line.split(' ')[0])
        velocity[1] = float(line.split(' ')[1])
        velocity[2] = float(line.split(' ')[2])
        node.SetValue(KratosMultiphysics.VELOCITY, velocity)


metric_parameters = KratosMultiphysics.Parameters("""
    {
        "minimal_size"                        : 0.00000000001,
        "maximal_size"                        : 10000000.0,
        "enforce_current"                     : false,
        "hessian_strategy_parameters":
        {
            "non_historical_metric_variable"  : true,
            "estimate_interpolation_error"    : false,
            "interpolation_error"             : 0.01
        },
        "anisotropy_remeshing"                : true
    }    """)

find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(main_model_part)
find_nodal_h.Execute()
hessian_metric = KratosMeshing.ComputeHessianSolMetricProcess(main_model_part,KratosMultiphysics.VELOCITY_X, metric_parameters)
hessian_metric.Execute()
hessian_metric = KratosMeshing.ComputeHessianSolMetricProcess(main_model_part,KratosMultiphysics.VELOCITY_Y, metric_parameters)
hessian_metric.Execute()

CreateGidControlOutput("metric_2d", main_model_part)

remesh_parameters = KratosMultiphysics.Parameters()

MmgProcess = KratosMeshing.MmgProcess2D(main_model_part,remesh_parameters)
MmgProcess.Execute()

find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(main_model_part)
find_nodal_h.Execute()

interp_error = metric_parameters["hessian_strategy_parameters"]["interpolation_error"].GetDouble()
CreateGidControlOutput("CAARC2d"+"_interperror"+str(interp_error), main_model_part)
KratosMultiphysics.ModelPartIO("CAARC2d"+"_interperror"+str(interp_error), KratosMultiphysics.IO.WRITE | KratosMultiphysics.IO.MESH_ONLY).WriteModelPart(main_model_part)
