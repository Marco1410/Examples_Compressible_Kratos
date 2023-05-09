import os
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
        mu.sort(reverse=True)
    return mu

def plot_mu_values(mu,name,xl,xu,yl,yu):
    mu_a = np.zeros(len(mu))
    mu_m = np.zeros(len(mu))
    fig  = plt.figure()
    fig.subplots_adjust(top=0.8)
    fig.set_figwidth(15.0)
    fig.set_figheight(10.0)
    for i in range(len(mu)):
        mu_a[i] = mu[i][0] * 180 / np.pi
        mu_m[i] = mu[i][1]
        fig = plt.plot(mu_m[i], mu_a[i], "o", markersize = 7.0, label = str(np.round(mu_a[i],1))+"º_"+str(np.round(mu_m[i],2)))
    fig = plt.title('Mu Values')
    fig = plt.axis([xl,xu+0.1,yl,yu])
    fig = plt.ylabel('Alpha')
    fig = plt.xlabel('Mach')
    fig = plt.grid()
    fig = plt.legend()
    fig = plt.tight_layout()
    fig = plt.savefig(name)


def save_mu_parameters(mu):
    import pickle
    archivo = open('mu_values.dat', 'wb')
    pickle.dump(mu, archivo)
    archivo.close()

if __name__ == "__main__":
    
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Data')
    os.mkdir("Data")
    
    with open("ProjectParameters_transonic.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    ################################
    ################################
    ################################
    def M_cr(critical_cp,critical_mach):
        f = critical_cp - ((1/(0.7*critical_mach**2))*(((2+0.4*critical_mach**2)/2.4)**3.5-1))
        return f

    def shock_parameters(mach_infinity,angle_of_attack):
        critical_cp = 1 - (1/mach_infinity**2)
        # Bisection method
        a = 0.2
        b = 1.0
        tol=1.0e-6
        critical_mach = (a + b) / 2.0
        while True:
            if b - a < tol:
                return x
            elif np.sign(M_cr(critical_cp,critical_mach)) * np.sign(M_cr(critical_cp,critical_mach)) > 0:
                a = critical_mach
            else:
                b = critical_mach
            critical_mach = (a + b) / 2.0

            if angle_of_attack > 4:
                critical_mach = critical_mach * 0.10
            elif angle_of_attack > 3:
                critical_mach = critical_mach * 0.20
            elif angle_of_attack > 2:
                critical_mach = critical_mach * 0.50
            elif angle_of_attack > 1:
                critical_mach = critical_mach * 0.70
            elif angle_of_attack > 0.5:
                critical_mach = critical_mach * 0.90

            upwind_factor_constant = 0.8
            return critical_mach,upwind_factor_constant
        
    mach_range  = [ 0.7,0.8]
    angle_range = [1.5,2.5]
    number_of_point_test = 25
    ################################
    ################################
    ################################

    mu_test  = get_multiple_params_by_Halton(number_of_point_test,angle_range,mach_range)
    save_mu_parameters(mu_test)
    plot_mu_values(mu_test,"Mu_Values.png",mach_range[0],mach_range[1],angle_range[0],angle_range[1])

    mesh_list = ["mesh_1"]#,"mesh_2"]

    # airfoils vs cp vs x
    fig  = plt.figure()
    fig.subplots_adjust(top=0.8)
    fig.set_figwidth(20.0)
    fig.set_figheight(15.0)

    ax1 = fig.add_subplot(211)
    ax1.axis([-0.1,1.15,-3.0,1.5])
    ax1.set_ylabel('Cp')
    ax1.set_title('Cp vs x')
    ax1.invert_yaxis()
    ax1.grid()

    ax2 = fig.add_subplot(212)
    ax2.axis([-0.1,1.1,-0.2,0.2])
    ax2.set_ylabel('y')
    ax2.set_xlabel('x')
    ax2.set_title('Airfoils')
    ax2.grid()

    for mesh_name in mesh_list:
        parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(mesh_name+"_Kratos")
        for id, mu in enumerate(mu_test):
            angle_of_attack = mu[0]
            mach_infinity = mu[1]
            critical_mach,upwind_factor_constant = shock_parameters(mach_infinity,angle_of_attack)
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
            ax1.plot( x, cp, "o", markersize = 7.0, label = mesh_name + "_" + str(np.round(mu[0] * 180 / np.pi,1))+"º_" + str(np.round(mu[1],2)))
        ax2.plot( x, y, "o", markersize = 5.0, label = mesh_name)

    ax1.legend()
    ax2.legend()
    fig = plt.tight_layout()
    fig = plt.savefig("Airfoils_Cp_x.png")
    

