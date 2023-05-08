from scipy.stats import qmc
import math
import numpy as np
import matplotlib.pyplot as plt
import KratosMultiphysics
from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis import PotentialFlowAnalysis

def get_multiple_params_by_Halton(number_of_values):
    sampler = qmc.Halton(d=2)
    sample = sampler.random(number_of_values)
    #Angle of attack
    l_angle = 0.0
    u_angle = 5.0
    #Mach infinit
    l_mach = 0.7
    u_mach = 0.9
    mu = []
    values = qmc.scale(sample, [l_angle,l_mach], [u_angle,u_mach])
    for i in range(number_of_values):
        #Angle of attack , Mach infinit
        mu.append([values[i,0] * math.pi / 180.0, values[i,1]])
    return mu

def plot_mu_values(mu,name):
    mu_a = np.zeros(len(mu))
    mu_m = np.zeros(len(mu))
    for i in range(len(mu)):
        mu_a[i] = mu[i][0] * 180 / math.pi
        mu_m[i] = mu[i][1]
    plt.plot(mu_m, mu_a, 'bs')
    plt.title('Mu Values')
    plt.ylabel('Alpha')
    plt.xlabel('Mach')
    plt.grid()
    plt.legend()
    plt.legend(bbox_to_anchor=(.85, 1.03, 1., .102), loc='upper left', borderaxespad=0.)
    plt.savefig(name)

if __name__ == "__main__":
    
    with open("ProjectParameters_transonic.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    mu_test  = get_multiple_params_by_Halton(4)
    plot_mu_values(mu_test,"Mu_Values.png")

    mesh = ["mesh_1_Kratos","mesh_2_Kratos"]
    for mesh_name in mesh:
        parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(mesh_name)
        for id, mu in enumerate(mu_test):

            parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(mu[0])
            parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mu[1])

            
            
            parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["critical_mach"].SetDouble(mu[1])
            parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["upwind_factor_constant"].SetDouble(mu[1])

            model = KratosMultiphysics.Model()
            simulation = PotentialFlowAnalysis(model,parameters)
            simulation.Run()

            modelpart = model["FluidModelPart.walls"]
            fout=open("Data/transonic_walls_"+mesh_name+"_"+str(id)+".dat",'w')
            fout.write("#  CoordinatesX CoordinatesY CoordinatesZ CoefPressure\n")
            for node in modelpart.Nodes:
                x=node.X0 ; y=node.Y0 ; z=node.Z0
                val=node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
                fout.write("%s %s %s %s\n" %(x,y,z,val))
            fout.close()
    

