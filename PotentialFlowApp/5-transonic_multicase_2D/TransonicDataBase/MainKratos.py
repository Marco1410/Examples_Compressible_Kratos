import time 
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis import PotentialFlowAnalysis


class PotentialFlowAnalysisWithFlush(PotentialFlowAnalysis):

    def __init__(self,model,project_parameters):
        super().__init__(model,project_parameters)

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def Clean():
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Results')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Data')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Captures')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('TransonicDataBase.post.lst')
    os.mkdir("Results")
    os.mkdir("Data")
    os.mkdir("Captures")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


if __name__ == "__main__":
    
    Clean()

    time_0 = time.time()

    with open("ProjectParameters_transonic.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    simulation_number = 0

    angles = np.arange( 1.50, 2.001, 0.02)
    machs  = np.arange( 0.72, 0.761, 0.005)

    critical_mach          = 0.9
    upwind_factor_constant = 1.5

    snapshot_variables_list = [CPFApp.VELOCITY_POTENTIAL,CPFApp.AUXILIARY_VELOCITY_POTENTIAL]

    for angle_of_attack in angles:

        angle_of_attack = np.round(angle_of_attack,2)

        for mach_infinity in machs:
            
            cp_min = 0
            cp_max = 0
            # cp vs x
            fig = plt.figure()
            fig.set_figwidth(12.0)
            fig.set_figheight(8.0)
            plt.ylabel('Cp')
            plt.xlabel('x')
            plt.title('Cp vs x - ' + " angle: " + str(angle_of_attack) + "º")
            plt.grid()

            mach_infinity = np.round(mach_infinity,3)

            parameters["solver_settings"]["scheme_settings"]["critical_mach"].SetDouble(critical_mach)
            parameters["solver_settings"]["scheme_settings"]["upwind_factor_constant"].SetDouble(upwind_factor_constant)
            parameters["solver_settings"]["scheme_settings"]["update_critical_mach"].SetDouble(critical_mach)
            parameters["solver_settings"]["scheme_settings"]["update_upwind_factor_constant"].SetDouble(upwind_factor_constant)
            parameters["solver_settings"]["scheme_settings"]["update_transonic_tolerance"].SetDouble(1e-30)
            parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach_infinity)
            parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(angle_of_attack)
            parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString("Results/" + "FOM_" + str(angle_of_attack) + "_" + str(mach_infinity))
            model = KratosMultiphysics.Model()
            simulation = PotentialFlowAnalysisWithFlush(model,parameters)
            simulation.Run()

            modelpart = model["FluidModelPart.Body2D_Body"]
            fout=open("Data/" + "FOM_" + str(angle_of_attack) + "_" + str(mach_infinity) + ".dat",'w')
            x    = np.zeros(modelpart.NumberOfNodes())
            y    = np.zeros(modelpart.NumberOfNodes())
            z    = np.zeros(modelpart.NumberOfNodes())
            cp   = np.zeros(modelpart.NumberOfNodes())
            for i,node in enumerate(modelpart.Nodes):
                x[i] = node.X0
                y[i] = node.X0
                z[i] = node.X0
                cp[i] = node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
                fout.write("%s %s %s %s\n" %(x[i],y[i],z[i],cp[i]))
            fout.close()
            plt.plot( x, cp, "xr", markersize = 3.5, label = "mach: " + str(mach_infinity))
            if np.min(cp) < cp_min:
                cp_min = np.min(cp)
            if np.max(cp) > cp_max:
                cp_max = np.max(cp)
            plt.axis([-0.05,1.35,cp_max+0.1,cp_min-0.1])
            plt.legend()
            plt.tight_layout()
            plt.savefig("Captures/" + str(simulation_number) + "_A_" + str(angle_of_attack) + "_M_" + str(mach_infinity) + ".png")
            plt.close(fig)
            
            simulation_number += 1

    time_final = time.time()

    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print(":::::: USED TIME:", np.round(time_final - time_0,2))
    print(":::::: SIMULATIONS NUMBER:", simulation_number)
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
