import sys
import time
import importlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

import KratosMultiphysics

def CreateAnalysisStageWithFlushInstance(cls, global_model, parameters):
    class AnalysisStageWithFlush(cls):

        def __init__(self, model,project_parameters, flush_frequency=10.0):
            super().__init__(model,project_parameters)
            self.flush_frequency = flush_frequency
            self.last_flush = time.time()
            sys.stdout.flush()

        def Initialize(self):
            super().Initialize()
            sys.stdout.flush()

        def FinalizeSolutionStep(self):
            super().FinalizeSolutionStep()

            if self.parallel_type == "OpenMP":
                now = time.time()
                if now - self.last_flush > self.flush_frequency:
                    sys.stdout.flush()
                    self.last_flush = now

    return AnalysisStageWithFlush(global_model, parameters)

if __name__ == "__main__":

    with open("ProjectParameters.json", 'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    analysis_stage_module_name = parameters["analysis_stage"].GetString()
    analysis_stage_class_name = analysis_stage_module_name.split('.')[-1]
    analysis_stage_class_name = ''.join(x.title() for x in analysis_stage_class_name.split('_'))

    analysis_stage_module = importlib.import_module(analysis_stage_module_name)
    analysis_stage_class = getattr(analysis_stage_module, analysis_stage_class_name)


    mach_infinity = 0.7
    angle_of_attack = 4 * np.pi / 180
    upwind_factor_constant = 2.5
    critical_mach = 0.6
    mach_number_limit = 1.45
    
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach_infinity)
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["critical_mach"].SetDouble(critical_mach)
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_number_limit"].SetDouble(mach_number_limit)
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(angle_of_attack)
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["upwind_factor_constant"].SetDouble(upwind_factor_constant)
    
    global_model = KratosMultiphysics.Model()
    simulation = CreateAnalysisStageWithFlushInstance(analysis_stage_class, global_model, parameters)
    simulation.Run()

    modelpart = global_model["FluidModelPart.Body2D_Body"]
    mainmodelpart = global_model["FluidModelPart"]
    KratosMultiphysics.ModelPartIO("new_meshe", KratosMultiphysics.IO.WRITE | KratosMultiphysics.IO.MESH_ONLY).WriteModelPart(mainmodelpart)

    x = np.zeros(modelpart.NumberOfNodes())
    y = np.zeros(modelpart.NumberOfNodes())
    z = np.zeros(modelpart.NumberOfNodes())
    cp = np.zeros(modelpart.NumberOfNodes())
    rho = np.zeros(modelpart.NumberOfNodes())
    for i,node in enumerate(modelpart.Nodes):
        x[i] = node.X0 ; y[i] = node.Y0 ; z[i] = node.Z0
        cp[i] = node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
    
    # Plot cp vs x
    fig,ax  = plt.subplots()
    fig.set_figwidth(15.0)
    fig.set_figheight(10.0)
    ax.plot( x, cp, "o", markersize = 3.0)
    ax.grid()
    plt.ylabel('Cp')
    plt.xlabel('x')
    plt.title('Cp vs x')
    ax.invert_yaxis()
    plt.tight_layout()
    fig.savefig("Airfoils_Cp_x.png")