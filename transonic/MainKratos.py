from scipy.stats import qmc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import KratosMultiphysics
from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis import PotentialFlowAnalysis

def get_multiple_params_by_Halton(number_of_values,angle,mach):
    sampler = qmc.Halton(d=2)
    sample = sampler.random(number_of_values)
    mu = []
    values = qmc.scale(sample, [angle[0],mach[0]], [angle[1],mach[1]])
    for i in range(number_of_values):
        #Angle of attack , Mach infinit
        mu.append([values[i,0] * np.pi / 180.0, values[i,1]])
    return mu

def plot_mu_values(mu,name):
    mu_a = np.zeros(len(mu))
    mu_m = np.zeros(len(mu))
    for i in range(len(mu)):
        mu_a[i] = mu[i][0] * 180 / np.pi
        mu_m[i] = mu[i][1]
        plt.plot(mu_m[i], mu_a[i], "x", label = str(np.round(mu_m[i],2))+"_"+str(np.round(mu_a[i],1))+"º")
    plt.title('Mu Values')
    plt.ylabel('Alpha')
    plt.xlabel('Mach')
    plt.grid()
    plt.legend()
    plt.legend()
    plt.savefig(name)

if __name__ == "__main__":
    
    with open("ProjectParameters_transonic.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    mu_test  = get_multiple_params_by_Halton(10,[0.0,1.0],[0.7,0.9])
    plot_mu_values(mu_test,"Mu_Values.png")

    mesh_list = ["mesh_1","mesh_2"]

    # airfoils ...............
    fig1, ax1 = plt.subplots()
    fig1.set_figwidth(6.0)
    fig1.set_figheight(3.4)

    # cps ...............
    fig2, ax2 = plt.subplots()
    fig2.set_figwidth(12.0)
    fig2.set_figheight(6.8)

    for mesh_name in mesh_list:
        parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(mesh_name+"_Kratos")
        for id, mu in enumerate(mu_test):
            angle_of_attack = mu[0]
            mach_infinity = mu[1]
            critical_mach = 0.9
            upwind_factor_constant = 1.0
            parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(angle_of_attack)
            parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach_infinity)
            parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["critical_mach"].SetDouble(critical_mach)
            parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["upwind_factor_constant"].SetDouble(upwind_factor_constant)

            model = KratosMultiphysics.Model()
            simulation = PotentialFlowAnalysis(model,parameters)
            simulation.Run()

            modelpart = model["FluidModelPart.walls"]
            fout=open("Data/"+mesh_name+"_"+str(id)+".dat",'w')
            fout.write("#  CoordinatesX CoordinatesY CoordinatesZ CoefPressure\n")
            x = np.zeros(modelpart.NumberOfNodes())
            y = np.zeros(modelpart.NumberOfNodes())
            z = np.zeros(modelpart.NumberOfNodes())
            cp = np.zeros(modelpart.NumberOfNodes())
            for i,node in enumerate(modelpart.Nodes):
                x[i] = node.X0 ; y[i] = node.Y0 ; z[i] = node.Z0
                cp[i] = node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
                fout.write("%s %s %s %s\n" %(x[i],y[i],z[i],cp[i]))
            fout.close()
            ax2.plot( x, cp, "x", markersize = 3.0, label = mesh_name + "_" + str(np.round(mu[1],2))+"_"+str(np.round(mu[0] * 180 / np.pi,1))+"º")
        ax1.plot( x, y, "+", markersize = 3.0, label = mesh_name)

    ax1.grid()
    plt.title("Airfoils")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.tight_layout()
    plt.legend()
    fig1.savefig("Airfoils.png")
    
    ax2.grid()
    plt.title("Cp vs x")
    plt.xlabel("x")
    plt.ylabel("cp")
    ax2.invert_yaxis()
    plt.tight_layout()
    plt.legend()
    fig2.savefig("Cp_vs_x.png")
    

