import sys
import time
import importlib
import numpy as np
import matplotlib.pyplot as plt
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
import matplotlib
matplotlib.use('Agg')
import KratosMultiphysics

def M_cr(cp,critical_mach):
    return cp-((1/(0.7*critical_mach**2))*((((2+0.4*critical_mach**2)/2.4)**3.5)-1))

def transonic_parameters(mach_infinity):
    # -0.6 cp min del naca 0012 a 0º
    # debe ser el menor cp del ala/perfil en regimen incompresible
    cp_min = -0.6

    # Karman-Tsien rule
    cp_KT = cp_min/(np.sqrt(1-mach_infinity**2)+(cp_min*mach_infinity**2/(1+np.sqrt(1-mach_infinity**2))))

    # Bisection method
    a = 0.02
    b = 0.95
    tol=1.0e-6
    critical_mach = (a + b) / 2.0
    while True:
        if b - a < tol:
            return np.round(critical_mach,2)
        elif M_cr(cp_KT,a) * M_cr(cp_KT,critical_mach) > 0:
            a = critical_mach
        else:
            b = critical_mach
        critical_mach = (a + b) / 2.0

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
        
        def InitializeSolutionStep(self):
            super().InitializeSolutionStep()

        def FinalizeSolutionStep(self):
            super().FinalizeSolutionStep()

            modelpart = global_model["FluidModelPart"]

            parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name" : "FluidModelPart",
                "shock_sensor"    : true,
                "artificial_bulk_viscosity_constant": 1.0,
                "calculate_nodal_area_at_each_step" : false
            }  """ )
            
            physics = CPFApp.ShockCapturingPhysicsProcess(modelpart, parameters)
            physics.Execute()

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


    mach_infinity   = 0.8
    angle_of_attack = -0.01
    critical_mach   = transonic_parameters(mach_infinity)
    mach_number_limit      = np.sqrt(3)

    upwind_factor_constant = 0.8

    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach_infinity)
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["critical_mach"].SetDouble(critical_mach)
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_number_limit"].SetDouble(mach_number_limit)
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(angle_of_attack * np.pi / 180)
    
    convergence_ratio = 1.0
    step = 0
    reset = False

    while (convergence_ratio > 1e-10):

        if reset:
            parameters["solver_settings"]["maximum_iterations"].SetInt(500)
        else:
            parameters["solver_settings"]["maximum_iterations"].SetInt(80)

        print("::::::::::::::::::::::::::::::::::::::::::::::::::::::NEW SIMULATION ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
        parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["upwind_factor_constant"].SetDouble(upwind_factor_constant)
        global_model = KratosMultiphysics.Model()
        simulation = CreateAnalysisStageWithFlushInstance(analysis_stage_class, global_model, parameters)
        simulation.Run()

        convergence_ratio = global_model["FluidModelPart"].ProcessInfo[KratosMultiphysics.CONVERGENCE_RATIO]

        if convergence_ratio > 0.5:
            reset = False
            step += 1
            upwind_factor_constant += 0.1
            print("New UPWIND:", upwind_factor_constant, "  ::::::::::::::::::::::::::::::::::::::::")
            print("Convergence ratio:", convergence_ratio, "  ::::::::::::::::::::::::::::::::::::::::")
        else:
            reset = True
            print(":::::::::::::::::::::::::::::: RESTART ::::::::::::::::::::::::::::::::")

        if upwind_factor_constant > 5.0:
            KratosMultiphysics.Logger.Print("               NON CONVERGENCE")
            sys.exit('               STOP')

    print("FINAL UPWIND:", upwind_factor_constant, "  ::::::::::::::::::::::::::::::::::::::::")
    print("FINAL STEPS :", step, "  ::::::::::::::::::::::::::::::::::::::::")

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