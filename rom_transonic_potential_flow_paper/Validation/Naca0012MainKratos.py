#####################################################################################################
#####################################################################################################
######################                NACA VALIDATION
#####################################################################################################
######################  This data is from:
######################  Volpe, G., and Jameson, A., “Transonic Potential Flow Calculations
######################  by Two Articial Density Methods,” AIAA Journal, Vol. 26, No. 4,
######################  April 1988, pp. 425–429.
######################  doi:10.2514/3.9910
#####################################################################################################

import KratosMultiphysics
import sys
import os
import time
import importlib
import numpy as np
import matplotlib.pyplot as plt
import KratosMultiphysics.FluidDynamicsApplication 

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

    mu_validation = []
    mu_validation.append([1.0,0.72])
    mu_validation.append([1.0,0.73])
    mu_validation.append([1.0,0.75])
    mu_validation.append([2.0,0.75])

    with open("Naca0012ProjectParameters.json", 'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    analysis_stage_module_name = parameters["analysis_stage"].GetString()
    analysis_stage_class_name = analysis_stage_module_name.split('.')[-1]
    analysis_stage_class_name = ''.join(x.title() for x in analysis_stage_class_name.split('_'))

    analysis_stage_module = importlib.import_module(analysis_stage_module_name)
    analysis_stage_class = getattr(analysis_stage_module, analysis_stage_class_name)

    fig = plt.figure(figsize=(10, 10))

    positions = [
    [0.075, 0.560, 0.4, 0.4],
    [0.575, 0.560, 0.4, 0.4],
    [0.075, 0.075, 0.4, 0.4],
    [0.575, 0.075, 0.4, 0.4]
    ]

    for n, mu in enumerate(mu_validation):

        angle_of_attack = mu[0]
        mach_infinity   = mu[1]

        case_name = f'{angle_of_attack}, {mach_infinity}'

        parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(np.double(angle_of_attack))
        parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(np.double(mach_infinity))
        parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString(f'Results_Naca0012/{angle_of_attack}, {mach_infinity}')
        parameters["output_processes"]["vtk_output"][0]["Parameters"]["output_path"].SetString(f'Results_Naca0012/{angle_of_attack}, {mach_infinity}')

        if (mach_infinity > 0.73 or angle_of_attack > 1.00):
            parameters["solver_settings"]["scheme_settings"]["initial_critical_mach"].SetDouble(0.90)
            parameters["solver_settings"]["scheme_settings"]["initial_upwind_factor_constant"].SetDouble(3.0)
            parameters["solver_settings"]["scheme_settings"]["update_relative_residual_norm"].SetDouble(1e-3)
            parameters["solver_settings"]["scheme_settings"]["target_critical_mach"].SetDouble(0.95)
            parameters["solver_settings"]["scheme_settings"]["target_upwind_factor_constant"].SetDouble(2.5)

        if (mach_infinity > 0.74 or angle_of_attack > 1.75):
            parameters["solver_settings"]["scheme_settings"]["initial_critical_mach"].SetDouble(0.80)
            parameters["solver_settings"]["scheme_settings"]["initial_upwind_factor_constant"].SetDouble(4.0)
            parameters["solver_settings"]["scheme_settings"]["update_relative_residual_norm"].SetDouble(1e-3)
            parameters["solver_settings"]["scheme_settings"]["target_critical_mach"].SetDouble(0.90)
            parameters["solver_settings"]["scheme_settings"]["target_upwind_factor_constant"].SetDouble(2.0)

        if (mach_infinity > 0.75 or angle_of_attack > 2.00):
            parameters["solver_settings"]["scheme_settings"]["initial_critical_mach"].SetDouble(0.80)
            parameters["solver_settings"]["scheme_settings"]["initial_upwind_factor_constant"].SetDouble(8.0)
            parameters["solver_settings"]["scheme_settings"]["update_relative_residual_norm"].SetDouble(1e-3)
            parameters["solver_settings"]["scheme_settings"]["target_critical_mach"].SetDouble(0.85)
            parameters["solver_settings"]["scheme_settings"]["target_upwind_factor_constant"].SetDouble(7.0)

        global_model = KratosMultiphysics.Model()
        simulation = CreateAnalysisStageWithFlushInstance(analysis_stage_class, global_model, parameters)
        simulation.Run()

        sub_ax = fig.add_axes(positions[n])

        #### Validation data ##########################################################
        validation_skin_data_filename = f"../reference_data/flo36/{case_name}.dat"
        if os.path.exists(validation_skin_data_filename):
            x_validation  = np.loadtxt(validation_skin_data_filename, usecols=(0,))
            cp_validation = np.loadtxt(validation_skin_data_filename, usecols=(1,))
            sub_ax.plot(x_validation, -cp_validation, 'r.-', markersize = 5.0, label = 'flo36')
        ###############################################################################

        modelpart = global_model["MainModelPart.Body2D_Body"]
        x = np.zeros(modelpart.NumberOfNodes())
        y = np.zeros(modelpart.NumberOfNodes())
        cp = np.zeros(modelpart.NumberOfNodes())
        for i,node in enumerate(modelpart.Nodes):
            x[i] = node.X; y[i] = node.Y
            cp[i] = node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
        y_gtz = y > 0
        y_lez = y <= 0
        x_gtz  = x[y_gtz] 
        cp_gtz = cp[y_gtz]
        indices_ordenados = sorted(range(len(x_gtz)), key=lambda i: x_gtz[i])
        x_gtz = [x_gtz[i] for i in indices_ordenados]
        cp_gtz = [cp_gtz[i] for i in indices_ordenados]
        x_lez  = x[y_lez] 
        cp_lez = cp[y_lez]
        indices_ordenados = sorted(range(len(x_lez)), key=lambda i: x_lez[i])
        x_lez = [x_lez[i] for i in indices_ordenados]
        cp_lez = [cp_lez[i] for i in indices_ordenados]
        cp_val_gtz = np.interp(x_gtz, x_validation[79:][::-1], cp_validation[79:][::-1])
        cp_val_lez = np.interp(x_lez, x_validation[:78], cp_validation[:78])
        x_sim = np.concatenate((x_lez,x_gtz[::-1]))
        cp_sim = np.concatenate((cp_lez, cp_gtz[::-1]))
        cp_val = np.concatenate((cp_val_lez,cp_val_gtz[::-1]))
        error = (np.linalg.norm(cp_sim-cp_val)/np.linalg.norm(cp_val))
        sub_ax.text(0.28, -1.0, f'Difference: {error:.2E}', size=12, bbox=dict(boxstyle="round", facecolor='red', alpha=0.5))
        sub_ax.plot( x, -cp, '.', markersize = 3.0, label=f'Kratos')

        sub_ax.set_title(f'Angle={angle_of_attack}, Mach={mach_infinity}')
        sub_ax.set_xlabel('x')
        sub_ax.set_ylabel('-Cp')
        sub_ax.grid()
        sub_ax.legend(loc='upper right')
        
    plt.savefig("Naca0012 Validation.pdf", dpi=400)
    # plt.show()
    plt.close()