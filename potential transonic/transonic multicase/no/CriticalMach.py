import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import KratosMultiphysics
from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis import PotentialFlowAnalysis

def M_cr(cp_min,critical_mach):
    #Karman-Tsien rule (mejor que Prandtl-Glauert que no es preciso para M > 0.7)    +(cp_min*critical_mach**2/(1+np.sqrt(1-critical_mach**2)))
    return cp_min/(np.sqrt(1-critical_mach**2))-((1/(0.7*critical_mach**2))*((((2+0.4*critical_mach**2)/2.4)**3.5)-1))

def shock_parameters(mach_infinity,angle_of_attack):
    # -0.6 cp min del naca 0012 a 0º
    # debe ser el menor cp del ala/perfil en regimen incompresible
    cp_min = -0.6
    # Bisection method
    a = 0.02
    b = 0.95
    tol=1.0e-6
    critical_mach = (a + b) / 2.0
    while True:
        if b - a < tol:
            critical_mach = 0.5
            mach_number_limit = 1.7320
            #upwind_factor_constant = 2.85
            return np.round(critical_mach,2),np.round(mach_number_limit,2)
        elif M_cr(cp_min,a) * M_cr(cp_min,critical_mach) > 0:
            a = critical_mach
        else:
            b = critical_mach
        critical_mach = (a + b) / 2.0

if __name__ == "__main__":
    
    with open("ProjectParameters_transonic.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    mesh_list = ["naca0012_0aoa"]
    angles = np.arange(3.0,4.5,0.5)
    machs = np.arange(0.60,0.90,0.05)
    upwind_factor_constants = np.arange(1.0,5.5,0.5)

    count = 0
    for mesh_name in mesh_list:
    
        for angle_of_attack in angles:
            angle_of_attack = np.round(angle_of_attack,2)

            for mach_infinity in machs:
                mach_infinity = np.round(mach_infinity,2)

                # #plot cp en 3D
                cp3d = plt.figure()
                cp3d.subplots_adjust()
                cp3d.set_figwidth(50.0)
                cp3d.set_figheight(25.0)
                # Creamos el plano 3D
                ax3 = cp3d.add_subplot(111, projection='3d')
                ax3.set_xlabel('x')
                ax3.set_ylabel('y')
                ax3.set_zlabel('Cp')    
                ax3.set_title('Cp vs x - ' + mesh_name + " - angle: " + str(angle_of_attack) + "º mach: " + str(mach_infinity))
                ax3.invert_zaxis()
                ax3.grid()
        
                for id, upwind in enumerate(upwind_factor_constants):
                    upwind = np.round(upwind,2)
                
                    critical_mach,mach_number_limit = shock_parameters(mach_infinity,angle_of_attack)

                    parameters["output_processes"].RemoveValue("gid_output")
                    parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(mesh_name)
                    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(angle_of_attack*np.pi/180)
                    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach_infinity)
                    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["critical_mach"].SetDouble(critical_mach)
                    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_number_limit"].SetDouble(mach_number_limit)
                    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["upwind_factor_constant"].SetDouble(upwind)

                    model = KratosMultiphysics.Model()
                    simulation = PotentialFlowAnalysis(model,parameters)
                    simulation.Run()

                    modelpart = model["FluidModelPart.Body2D_Body"]
                    x = np.zeros(modelpart.NumberOfNodes())
                    cp = np.zeros(modelpart.NumberOfNodes())
                    for i,node in enumerate(modelpart.Nodes):
                        x[i] = node.X0
                        cp[i] = node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)

                    ax3.scatter(x, id*10, cp, label = "critical mach: " + str(critical_mach) + " mach limit: " + str(mach_number_limit) + " Upwind: " + str(upwind))
                
                ax3.legend()
                cp3d = plt.tight_layout()
                cp3d = plt.savefig("plots/Airfoils Cp vs x - " + mesh_name + " - " + " angle: " + str(angle_of_attack) + "º mach: " + str(mach_infinity) + "_" + str(count) + ".png")
                count += 1