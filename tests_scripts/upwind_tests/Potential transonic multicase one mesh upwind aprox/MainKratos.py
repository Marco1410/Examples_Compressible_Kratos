import time 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
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

def upwind_aprox(angle,mach): #0.0 0.6 / 5.0 0.85
    return 2.1489732 /(1 + np.exp(-(1.55470776*angle+37.876707*mach-31.48703546))) + 0.09767296 + 0.5

if __name__ == "__main__":
    
    time_0 = time.time()

    with open("ProjectParameters_transonic.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    simulation_number = 0

    angles = np.arange( 0.00, 5.01 , 0.5)
    machs  = np.arange( 0.60, 0.801, 0.05)

    for angle_of_attack in angles:

        angle_of_attack = np.round(angle_of_attack,2)

        for mach_infinity in machs:
            
            #PLOT:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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

            #print("::::::::::::::::::::::::::::::::::::NEW CASE::::::::::::::::::::::::::::::::::::")
            simulation_number += 1

            mach_infinity = np.round(mach_infinity,2)
            critical_mach = transonic_parameters(mach_infinity)
            upwind_factor_constant = upwind_aprox(angle_of_attack,mach_infinity)
            parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(-0.025 + angle_of_attack*np.pi/180)
            parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach_infinity)
            parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["critical_mach"].SetDouble(critical_mach)
            parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["upwind_factor_constant"].SetDouble(np.round(upwind_factor_constant,3))
            
            model = KratosMultiphysics.Model()
            simulation = PotentialFlowAnalysisWithFlush(model,parameters)
            simulation.Run()

            modelpart = model["FluidModelPart.Body2D_Body"]
            x    = np.zeros(modelpart.NumberOfNodes())
            cp   = np.zeros(modelpart.NumberOfNodes())
            for i,node in enumerate(modelpart.Nodes):
                x[i] = node.X0
                cp[i] = node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
            plt.plot( x, cp, "xr", markersize = 3.5, label = "mach: " + str(mach_infinity) + " critical mach: " + str(np.round(critical_mach,2)) + " upwind: " + str(np.round(upwind_factor_constant,3)))
            if np.min(cp) < cp_min:
                cp_min = np.min(cp)
            if np.max(cp) > cp_max:
                cp_max = np.max(cp)
            plt.axis([-0.05,1.35,cp_max+0.1,cp_min-0.1])
            plt.legend()
            plt.tight_layout()
            plt.savefig("Simulation_" + str(simulation_number) + "_A_" + str(angle_of_attack) + "_M_" + str(mach_infinity) + ".png")
            plt.close(fig)

    time_final = time.time()
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print(":::::: USED TIME:", np.round(time_final - time_0,2))
    print(":::::: SIMULATIONS NUMBER:", simulation_number)
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
