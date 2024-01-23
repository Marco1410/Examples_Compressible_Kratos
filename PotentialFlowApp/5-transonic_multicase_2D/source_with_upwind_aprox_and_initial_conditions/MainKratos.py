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

    def ModifyAfterSolverInitialize(self):
        if (os.path.isfile("previous_data_1_2.dat")):
            id_nodes = np.loadtxt("previous_data_1_2.dat",usecols=(0,))
            velocity_potential = np.loadtxt("previous_data_1_2.dat",usecols=(1,))
            auxiliary_velocity_potential = np.loadtxt("previous_data_1_2.dat",usecols=(2,))
            for i in range(len(id_nodes)):
                node = model["FluidModelPart"].GetNode(int(id_nodes[i]))
                node.SetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL, 0, velocity_potential[i])
                node.SetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL, 1, velocity_potential[i])
                node.SetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL, 0, auxiliary_velocity_potential[i])
                node.SetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL, 1, auxiliary_velocity_potential[i])
        elif(os.path.isfile("previous_data_1.dat")):
            id_nodes = np.loadtxt("previous_data_1.dat",usecols=(0,))
            velocity_potential = np.loadtxt("previous_data_1.dat",usecols=(1,))
            auxiliary_velocity_potential = np.loadtxt("previous_data_1.dat",usecols=(2,))
            for i in range(len(id_nodes)):
                node = model["FluidModelPart"].GetNode(int(id_nodes[i]))
                node.SetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL, 0, velocity_potential[i])
                node.SetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL, 1, velocity_potential[i])
                node.SetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL, 0, auxiliary_velocity_potential[i])
                node.SetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL, 1, auxiliary_velocity_potential[i])

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

def upwind_aprox(angle,mach):
    # valores interpolados en el rango: 
    # angle: 0.0, 5.0
    # mach : 0.60, 0.85
    return 2.1489732 /(1 + np.exp(-(1.55470776*angle+37.876707*mach-31.48703546))) + 0.59767296

if __name__ == "__main__":
    
    time_0 = time.time()

    with open("ProjectParameters_transonic.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    simulation_number = 0

    angles = np.arange( 0.00, 3.01 , 0.1)
    machs  = np.arange( 0.60, 0.851, 0.01)

    for angle_of_attack in angles:

        KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('previous_data_1.dat')
        KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('previous_data_1_2.dat')
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
            upwind_factor_constant = upwind_aprox(angle_of_attack, mach_infinity)

            parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(angle_of_attack*np.pi/180)
            parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach_infinity)
            parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["critical_mach"].SetDouble(critical_mach)

            convergence_ratio = 1.0
            absolute_norm = 1.0
            tolerancia = 1e-10
            write_and_plot = True

            while (convergence_ratio > tolerancia and absolute_norm > tolerancia):
                
                parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["upwind_factor_constant"].SetDouble(np.round(upwind_factor_constant,3))
                model = KratosMultiphysics.Model()
                simulation = PotentialFlowAnalysisWithFlush(model,parameters)
                simulation.Run()
                
                convergence_ratio = model["FluidModelPart"].ProcessInfo[KratosMultiphysics.CONVERGENCE_RATIO]
                absolute_norm = model["FluidModelPart"].ProcessInfo[KratosMultiphysics.RESIDUAL_NORM]

                if convergence_ratio >= 0.5:
                    upwind_factor_constant += 0.2
                else:
                    upwind_factor_constant += 0.05
                    #print("::::::::::::::::::::::::::::::::NEW UPWIND::::::::::::::::::::::::::::::::::::::")
                
                if convergence_ratio < 0.5:
                    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('previous_data_1_2.dat')
                    fout=open("previous_data_1_2.dat",'w')
                    modelpart = model["FluidModelPart"]
                    for node in modelpart.Nodes:
                        id_node = node.Id
                        velocity_potential = node.GetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL,1)
                        auxiliary_velocity_potential = node.GetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL,1)
                        fout.write("%s %s %s\n" %(id_node, velocity_potential, auxiliary_velocity_potential))
                    fout.close()
                    #print(":::::::::::::::::::::::::::::::::RESTART / SAVE:::::::::::::::::::::::::::::::::")

                if upwind_factor_constant > 5.0:
                    print("::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
                    print("::::::::::::::::::::::::::::::::::::::: Non Convergence ::::::::::::::::::::::::::")
                    print("::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
                    write_and_plot = False
                    break

            KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('previous_data_1_2.dat')
            if write_and_plot:
                fout=open("previous_data_1.dat",'w')
                modelpart = model["FluidModelPart"]
                for node in modelpart.Nodes:
                    id_node = node.Id
                    velocity_potential = node.GetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL,1)
                    auxiliary_velocity_potential = node.GetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL,1)
                    fout.write("%s %s %s\n" %(id_node, velocity_potential, auxiliary_velocity_potential))
                fout.close()
            
                #print(":::::::::::::::::::::::::::::::::RESTART / SAVE:::::::::::::::::::::::::::::::::")

                #PLOT:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('previous_data_1.dat')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('previous_data_1_2.dat')
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print(":::::: USED TIME:", np.round(time_final - time_0,2))
    print(":::::: SIMULATIONS NUMBER:", simulation_number)
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
