import sys
import time
import importlib
import numpy as np
import matplotlib.pyplot as plt
import KratosMultiphysics.FluidDynamicsApplication 
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

            # modelpart = global_model["FluidModelPart"]
            # heat_capacity_ratio = modelpart.ProcessInfo.GetValue(KratosMultiphysics.FluidDynamicsApplication.HEAT_CAPACITY_RATIO)
            # for node in modelpart.Nodes:
            #     density  = node.GetValue(KratosMultiphysics.DENSITY)
            #     cp       = node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
            #     velocity = node.GetValue(KratosMultiphysics.VELOCITY)
            #     pressure = 0.5 * density * (velocity[1]**2+velocity[2]**2) * cp
            #     entropy  = (density / (heat_capacity_ratio - 1.0)) * np.log(np.abs(pressure) / density**heat_capacity_ratio)
            #     node.SetValue(KratosMultiphysics.FluidDynamicsApplication.NUMERICAL_ENTROPY,entropy)
    
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


    mach_infinity = 0.8
    angle_of_attack = 0.0 * np.pi / 180

    upwind_factor_constant = 2.5
    critical_mach = 0.7
    mach_number_limit = 1.25
    
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach_infinity)
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["critical_mach"].SetDouble(critical_mach)
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_number_limit"].SetDouble(mach_number_limit)
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(angle_of_attack)
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["upwind_factor_constant"].SetDouble(upwind_factor_constant)
    
    global_model = KratosMultiphysics.Model()
    simulation = CreateAnalysisStageWithFlushInstance(analysis_stage_class, global_model, parameters)
    simulation.Run()

    modelpart = global_model["FluidModelPart.Body2D_Body"]
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