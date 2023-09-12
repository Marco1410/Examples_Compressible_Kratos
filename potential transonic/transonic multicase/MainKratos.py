import time 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import KratosMultiphysics
from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis import PotentialFlowAnalysis

def M_cr(cp,critical_mach):
    return cp-((1/(0.7*critical_mach**2))*((((2+0.4*critical_mach**2)/2.4)**3.5)-1))

def transonic_parameters(mach_infinity,mesh_name):
    # -0.6 cp min del naca 0012 a 0º
    # debe ser el menor cp del ala/perfil en regimen incompresible
    # disminuye en valor absoluto con el espesor del perfil

    if mesh_name == "mesh_Kratos":
        cp_min = -0.55
    else:
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
    
    time_0 = time.clock_gettime(1)

    with open("ProjectParameters_transonic.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    mesh_list = ["mesh_Kratos", "naca0012_0aoa_00"]
    angles = np.arange(0.0,2.1,0.5)
    machs = np.arange(0.75,0.86,0.05)

    for angle_of_attack in angles:
        angle_of_attack = np.round(angle_of_attack,2)

        upwind_factor_constant = 0.1 
    
        for mesh_name in mesh_list:

            # # #plot cp en 3D
            # cp3d = plt.figure()
            # cp3d.subplots_adjust()
            # cp3d.set_figwidth(20.0)
            # cp3d.set_figheight(15.0)
            # # Creamos el plano 3D
            # ax3 = cp3d.add_subplot(111, projection='3d')
            # ax3.set_xlabel('x')
            # ax3.set_ylabel('y')
            # ax3.set_zlabel('Cp')    
            # ax3.set_title('Cp vs x - ' + mesh_name + " - angle: " + str(angle_of_attack) + "º" + "critical mach: " + str(critical_mach) + " mach limit: " + str(mach_number_limit) + " Upwind: " + str(upwind_factor_constant))
            # ax3.invert_zaxis()
            # ax3.grid()

            cp_min = 0
            cp_max = 0

            # airfoils vs cp vs x
            fig  = plt.figure()
            fig.subplots_adjust(top=0.8)
            fig.set_figwidth(20.0)
            fig.set_figheight(15.0)

            ax1 = fig.add_subplot(211)
            ax1.set_ylabel('Cp')
            ax1.set_title('Cp vs x - ' + mesh_name + " - angle: " + str(angle_of_attack) + "º")
            ax1.grid()

            ax2 = fig.add_subplot(212)
            ax2.set_ylabel('y')
            ax2.set_xlabel('x')
            ax2.axis([-0.05,1.35,-0.15,0.15])
            ax2.set_title('Airfoils')
            ax2.grid()

            for id,mach_infinity in enumerate(machs):
                mach_infinity = np.round(mach_infinity,2)
                critical_mach          = transonic_parameters(mach_infinity,mesh_name)
                mach_number_limit      = 1.7320 
                
                parameters["output_processes"].RemoveValue("gid_output")
                parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(mesh_name)
                parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(angle_of_attack*np.pi/180)
                parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach_infinity)
                parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["critical_mach"].SetDouble(critical_mach)
                parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_number_limit"].SetDouble(mach_number_limit)
    

                convergence_ratio = 1.0
                step = 0

                while (convergence_ratio > 1e-2):
                    print("NEW INTERNAL SIMULATION ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
                    print("MACH :", np.round(mach_infinity,2)," ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
                    print("ANGLE:", np.round(angle_of_attack,2)," ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
                    print("::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
                    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["upwind_factor_constant"].SetDouble(upwind_factor_constant)
        
                    model = KratosMultiphysics.Model()
                    simulation = PotentialFlowAnalysis(model,parameters)
                    simulation.Run()

                    convergence_ratio = model["FluidModelPart"].ProcessInfo[KratosMultiphysics.CONVERGENCE_RATIO]
                    if convergence_ratio > 1e-2:
                        step += 1
                        upwind_factor_constant += 0.1
                        print("New UPWIND:", upwind_factor_constant, "  ::::::::::::::::::::::::::::::::::::::::")
                        print("Convergence ratio:", convergence_ratio, "  ::::::::::::::::::::::::::::::::::::::::")

                        if upwind_factor_constant > 5.0:
                            KratosMultiphysics.Logger.Print("               NON CONVERGENCE")
                            # sys.exit('               STOP')
                            break


                modelpart = model["FluidModelPart.Body2D_Body"]
                x = np.zeros(modelpart.NumberOfNodes())
                y = np.zeros(modelpart.NumberOfNodes())
                cp = np.zeros(modelpart.NumberOfNodes())
                for i,node in enumerate(modelpart.Nodes):
                    x[i] = node.X0
                    y[i] = node.Y0
                    cp[i] = node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
                
                ax1.plot( x, cp, "o", markersize = 3.0, label = "mach: " + str(mach_infinity) + " critical mach: " + str(np.round(critical_mach,2)) + " mach limit: " + str(mach_number_limit) + " upwind: " + str(np.round(upwind_factor_constant,2)))

                if np.min(cp) < cp_min:
                    cp_min = np.min(cp)

                if np.max(cp) > cp_max:
                    cp_max = np.max(cp)
                # ax3.scatter(x, id*10, cp, label = "mach: " + str(mach_infinity))
                
            # ax3.legend()
            # cp3d = plt.tight_layout()
            # cp3d = plt.savefig("plots/Airfoils Cp vs x - " + mesh_name + " - " + " angle: " + str(angle_of_attack) + "º mach: " + str(mach_infinity) + "_" + str(count) + ".png")
            
            ax2.plot( x, y, "o", markersize = 3.0)
            
            ax1.axis([-0.05,1.35,cp_max+0.1,cp_min-0.1])
            ax1.legend()
            ax2.legend()
            fig = plt.tight_layout()
            fig = plt.savefig("plots/Airfoils Cp vs x " + mesh_name + " angle: " + str(angle_of_attack) + "º.png")
    time_final = time.clock_gettime(1)
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print(":::::: USED TIME:", time_final - time_0,"::::::::::::::::::")
    print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")