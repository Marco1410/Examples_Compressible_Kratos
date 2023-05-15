import os
import pickle
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
    values[0,1] = mach[0]
    values[number_of_values-1,1] = mach[1]
    for i in range(number_of_values):
        #Angle of attack , Mach infinit
        mu.append([values[i,1], values[i,0] * np.pi / 180.0])
        mu.sort(reverse=False)
    return mu

def plot_mu_values(mu,name,xl,xu,yl,yu):
    mu_a = np.zeros(len(mu))
    mu_m = np.zeros(len(mu))
    fig  = plt.figure()
    fig.subplots_adjust(top=0.8)
    fig.set_figwidth(20.0)
    fig.set_figheight(15.0)
    for i in range(len(mu)):
        mu_a[i] = mu[i][1] * 180 / np.pi
        mu_m[i] = mu[i][0]
        fig = plt.plot(mu_m[i], mu_a[i], "o", markersize = 5.0, label = str(np.round(mu_a[i],1))+"º_"+str(np.round(mu_m[i],2)))
    fig = plt.title('Mu Values')
    fig = plt.axis([xl,xu+0.1,yl,yu])
    fig = plt.ylabel('Alpha')
    fig = plt.xlabel('Mach')
    fig = plt.grid()
    fig = plt.legend()
    fig = plt.tight_layout()
    fig = plt.savefig(name)


def save_mu_parameters(mu):
    archivo = open('mu_values.dat', 'wb')
    pickle.dump(mu, archivo)
    archivo.close()

def load_mu_parameters():
    archivo = open('mu_values.dat', 'rb')
    lista = pickle.load(archivo)
    mu = np.asarray(lista)
    archivo.close()
    return mu

def M_cr(cp_min,critical_mach):
    #Karman-Tsien rule (mejor que Prandtl-Glauert que no es preciso para M > 0.7)
    return cp_min/np.sqrt(1-critical_mach**2)-((1/(0.7*critical_mach**2))*((((2+0.4*critical_mach**2)/2.4)**3.5)-1))

def shock_parameters(mach_infinity,angle_of_attack):
    # -0.6 cp min del naca 0012 a 0º
    # debe ser el menor cp del ala/perfil en regimen incompresible
    cp_min = -0.6
    # Bisection method
    a = 0.2
    b = 0.95
    tol=1.0e-6
    critical_mach = (a + b) / 2.0
    while True:
        if b - a < tol:
            #critical_mach = 0.6
            mach_number_limit = critical_mach * 2.0
            upwind_factor_constant = mach_number_limit**2 * 1.6
            print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
            print(mach_infinity,critical_mach,mach_number_limit,upwind_factor_constant)
            #input()
            return critical_mach,upwind_factor_constant,mach_number_limit
        elif M_cr(cp_min,a) * M_cr(cp_min,critical_mach) > 0:
            a = critical_mach
        else:
            b = critical_mach
        critical_mach = (a + b) / 2.0

if __name__ == "__main__":
    
    with open("ProjectParameters_transonic.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    ################################
    ################################
    ################################
    # en general el mach de crucero (en transonico) < 0.9
    # y el angulo de ataque < 5º
    mach_range  = [0.65,0.75]
    angle_range = [2.0,4.5]
    number_of_point_test = 15
    ################################
    ################################
    ################################

    #mu_test  = get_multiple_params_by_Halton(number_of_point_test,angle_range,mach_range)
    #save_mu_parameters(mu_test)
    #plot_mu_values(mu_test,"Mu_Values.png",mach_range[0],mach_range[1],angle_range[0],angle_range[1])
    mu_test = load_mu_parameters()

    mesh_list = ["naca0012_0aoa"]

    # #plot cp en 3D
    # cp3d = plt.figure()
    # cp3d.subplots_adjust()
    # cp3d.set_figwidth(15.0)
    # cp3d.set_figheight(10.0)
    # # Creamos el plano 3D
    # ax3 = cp3d.add_subplot(111, projection='3d')
    # ax3.set_xlabel('x')
    # ax3.set_ylabel('y')
    # ax3.set_zlabel('Cp')    
    # ax3.set_title('Cp vs x - naca0012')
    # ax3.invert_zaxis()
    # ax3.grid()

    # airfoils vs cp vs x
    fig  = plt.figure()
    fig.subplots_adjust(top=0.8)
    fig.set_figwidth(20.0)
    fig.set_figheight(15.0)

    ax1 = fig.add_subplot(211)
    ax1.axis([-0.6,0.65,-3.0,1.5])
    ax1.set_ylabel('Cp')
    ax1.set_title('Cp vs x')
    ax1.invert_yaxis()
    ax1.grid()

    ax2 = fig.add_subplot(212)
    ax2.axis([-0.6,0.65,-0.2,0.2])
    ax2.set_ylabel('y')
    ax2.set_xlabel('x')
    ax2.set_title('Airfoils')
    ax2.grid()

    for mesh_name in mesh_list:
        parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(mesh_name)
        fout=open("critical_cp_mach_"+mesh_name+".dat", 'w')
        fout.write("# Angle_of_atack mach_infinity Critical_Cp Critical_Mach upwind_factor_constant \n")
        for id, mu in enumerate(mu_test):
            angle_of_attack = mu[1]
            mach_infinity = mu[0]

            critical_mach,upwind_factor_constant,mach_number_limit = shock_parameters(mach_infinity,angle_of_attack)

            parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(angle_of_attack)
            parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach_infinity)
            parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["critical_mach"].SetDouble(critical_mach)
            parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["upwind_factor_constant"].SetDouble(upwind_factor_constant)
            parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_number_limit"].SetDouble(mach_number_limit)

            model = KratosMultiphysics.Model()
            simulation = PotentialFlowAnalysis(model,parameters)
            simulation.Run()

            modelpart = model["FluidModelPart.Body2D_Body"]
            x = np.zeros(modelpart.NumberOfNodes())
            y = np.zeros(modelpart.NumberOfNodes())
            z = np.zeros(modelpart.NumberOfNodes())
            cp = np.zeros(modelpart.NumberOfNodes())
            for i,node in enumerate(modelpart.Nodes):
                x[i] = node.X0 ; y[i] = node.Y0 ; z[i] = node.Z0
                cp[i] = node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
            ax1.plot( x, cp, "o", markersize = 3.0, label = mesh_name + "_" + str(np.round(mu[1] * 180 / np.pi,1))+"º_" + str(np.round(mu[0],2)))
            # ax3.scatter(x+0.5,id*10,cp,".",label=str(np.round(mu[1]*180/np.pi,1))+"º_"+str(np.round(mu[0],2)))
        fout.close()
        ax2.plot( x, y, "o", markersize = 3.0, label = mesh_name)

    ax1.legend()
    ax2.legend()
    fig = plt.tight_layout()
    fig = plt.savefig("Airfoils_Cp_x.png")
    # ax3.legend()
    # cp3d = plt.tight_layout()
    # cp3d = plt.savefig("cp_3d")
    

