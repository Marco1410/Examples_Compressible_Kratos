import time 
import numpy as np
import KratosMultiphysics
from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis import PotentialFlowAnalysis


class PotentialFlowAnalysisWithFlush(PotentialFlowAnalysis):

    def __init__(self,model,project_parameters):
        super().__init__(model,project_parameters)

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

def M_cr(cp,critical_mach):
    return cp-((1/(0.7*critical_mach**2))*((((2+0.4*critical_mach**2)/2.4)**3.5)-1))

def transonic_parameters(mach_infinity):
    # -0.6 cp min del naca 0012 a 0º
    # debe ser el menor cp del ala/perfil en regimen incompresible
    # disminuye en valor absoluto con el espesor del perfil
    cp_min = -0.6
    if mach_infinity < 0.7:
        # Prandtl-Glauert rule (no es preciso para M > 0.7)
        cp_M = cp_min/(np.sqrt(1-mach_infinity**2))
    else:
        # Karman-Tsien rule
        cp_M = cp_min/(np.sqrt(1-mach_infinity**2)+(cp_min*mach_infinity**2/(1+np.sqrt(1-mach_infinity**2))))
    # Bisection method
    a = 0.02
    b = 0.95
    tol=1.0e-6
    critical_mach = (a + b) / 2.0
    while True:
        if b - a < tol:
            return np.round(critical_mach,2)
        elif M_cr(cp_M,a) * M_cr(cp_M,critical_mach) > 0:
            a = critical_mach
        else:
            b = critical_mach
        critical_mach = (a + b) / 2.0

if __name__ == "__main__":
    
    time_0 = time.time()

    with open("ProjectParameters_transonic.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    simulation_number = 0

    angles = np.arange( 0.00, 3.01 , 0.5)
    machs  = np.arange( 0.60, 0.751, 0.05)

    for angle_of_attack in angles:

        angle_of_attack = np.round(angle_of_attack,2)

        for mach_infinity in machs:
            simulation_number += 1

            print("::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
            print(simulation_number,"::::::::::::::::::::::::::::::::::::::")
            print("::::::::::::::::::::::::::::::::::::::::::::::::::::::::")

            mach_infinity = np.round(mach_infinity,2)
            critical_mach = transonic_parameters(mach_infinity)
            parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString("GiD/simulation_" + str(simulation_number) + "_A_" + str(angle_of_attack) + "_M_"+str(mach_infinity))
            parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(angle_of_attack*np.pi/180)
            parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach_infinity)
            parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["critical_mach"].SetDouble(critical_mach)
            parameters["solver_settings"]["maximum_iterations"].SetInt(80)
            # parameters["output_processes"].RemoveValue("gid_output")

            upwind_factor_constant = 0.2
            convergence_ratio = 1.0
            absolute_norm = 1.0
            tolerancia = 1e-10

            while (convergence_ratio > tolerancia and absolute_norm > tolerancia):
                
                parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["upwind_factor_constant"].SetDouble(np.round(upwind_factor_constant,3))
                model = KratosMultiphysics.Model()
                simulation = PotentialFlowAnalysisWithFlush(model,parameters)
                simulation.Run()
                
                convergence_ratio = model["FluidModelPart"].ProcessInfo[KratosMultiphysics.CONVERGENCE_RATIO]
                absolute_norm = model["FluidModelPart"].ProcessInfo[KratosMultiphysics.RESIDUAL_NORM]

                if convergence_ratio >= 0.5:
                    upwind_factor_constant += 0.5
                else:
                    upwind_factor_constant += 0.05
                    parameters["solver_settings"]["maximum_iterations"].SetInt(300)
                
                if convergence_ratio < 1e-6:
                    break

                if upwind_factor_constant > 5.0:
                    print("::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
                    print("::::::::::::::::::::::::::::::: Non Convergence ::::::::::::::::::::::::::")
                    print("::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
                    break

    time_final = time.time()
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print(":::::: USED TIME:", np.round(time_final - time_0,2))
    print(":::::: SIMULATIONS NUMBER:", simulation_number)
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
