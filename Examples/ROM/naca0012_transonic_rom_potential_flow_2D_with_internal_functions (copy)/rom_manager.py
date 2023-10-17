import KratosMultiphysics
import numpy as np
import os
import pickle
from scipy.stats import qmc
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import KratosMultiphysics.kratos_utilities
from KratosMultiphysics.MeshMovingApplication.mesh_moving_analysis import MeshMovingAnalysis
from KratosMultiphysics.RomApplication.rom_manager import RomManager
    
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def GenerateMeshe(id, angle, typename):

    angle = -1.0 * angle * np.pi / 180.0 #Negativo porque se lo paso a la meshmoving app y lo paso a radianes
                
    with open("ProjectParametersMeshMoving.json",'r') as parameter_file:
        mesh_parameters = KratosMultiphysics.Parameters(parameter_file.read())
    mesh_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["rotation_angle"].SetDouble(angle)
    
    model = KratosMultiphysics.Model()
    mesh_simulation = MeshMovingAnalysis(model,mesh_parameters)
    mesh_simulation.Run()

    mainmodelpart = model["MainModelPart"]
    KratosMultiphysics.ModelPartIO("Meshes/" + typename + "_mesh_" + str(id), KratosMultiphysics.IO.WRITE | KratosMultiphysics.IO.MESH_ONLY).WriteModelPart(mainmodelpart)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def CustomizeSimulation(cls, global_model, parameters):

    class CustomSimulation(cls):

        def __init__(self, model,project_parameters, custom_param = None):
            super().__init__(model,project_parameters)
            self.custom_param  = custom_param

        def Initialize(self):
            super().Initialize()

        def FinalizeSolutionStep(self):
            super().FinalizeSolutionStep()

            if parameters["output_processes"].Has("gid_output"):
                nametype = parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].GetString()
                simulationtype = nametype.split('_')[0].split('/')[1]

                nametype = parameters["solver_settings"]["model_import_settings"]["input_filename"].GetString()
                typename = nametype.split('/')[-1].split('_')[0]
                id       = nametype.split('/')[1].split('_')[-1]
                # guardar aqui datos directamente de la skin y plotea
                if typename == "test":
                    fout=open("Data/"+simulationtype+"_"+typename+"_"+str(id)+".dat",'w')
                else:
                    fout=open("Data/"+simulationtype+"_"+typename+"_"+str(id)+".dat",'w')
                modelpart = self.model["MainModelPart.Body2D_Body"]
                for node in modelpart.Nodes:
                    x=node.X ; y=node.Y ; z=node.Z
                    cp=node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
                    fout.write("%s %s %s %s\n" %(x,y,z,cp))
                fout.close()

        def CustomMethod(self):
            return self.custom_param

    return CustomSimulation(global_model, parameters)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def UpdateProjectParameters(parameters, mu=None):
    mach_infinity          = mu[1]
    upwind_factor_constant = mu[2]
    id                     = mu[3]
    typename               = mu[4]
    critical_mach          = transonic_parameters(mach_infinity)

    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["upwind_factor_constant"].SetDouble(upwind_factor_constant)
    parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString("Meshes/" + typename + "_mesh_" + str(id))
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach_infinity)
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["critical_mach"].SetDouble(critical_mach)
    return parameters

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def UpdateMaterialParametersFile(material_parametrs_file_name, mu):
    pass

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def GetRomManagerParameters():
    general_rom_manager_parameters = KratosMultiphysics.Parameters("""{
            "rom_stages_to_train" : ["ROM"],      // ["ROM","HROM"]
            "rom_stages_to_test"  : ["ROM"],      // ["ROM","HROM"]
            "paralellism" : null,                        // null, TODO: add "compss"
            "projection_strategy": "lspg",               // "lspg", "galerkin", "petrov_galerkin"
            "assembling_strategy": "global",                                                       
            "save_gid_output": true,                     // false, true #if true, it must exits previously in the ProjectParameters.json
            "save_vtk_output": false,                    // false, true #if true, it must exits previously in the ProjectParameters.json
            "output_name": "mu",                         // "id" , "mu"
            "ROM":{
                "svd_truncation_tolerance": 1e-30,
                "model_part_name": "MainModelPart",                                      // This changes depending on the simulation: Structure, MainModelPart, ThermalPart #TODO: Idenfity it automatically
                "nodal_unknowns": ["VELOCITY_POTENTIAL","AUXILIARY_VELOCITY_POTENTIAL"], // Main unknowns. Snapshots are taken from these
                "rom_basis_output_format": "json",                                       // "json" "numpy"
                "rom_basis_output_name": "RomParameters",
                "snapshots_control_type": "step",                                        // "step", "time"
                "snapshots_interval": 1,
                "petrov_galerkin_training_parameters":{
                    "basis_strategy": "residuals",                                        // 'residuals', 'jacobian'
                    "include_phi": false,
                    "svd_truncation_tolerance": 1e-30,
                    "echo_level": 0
                },
                "lspg_rom_bns_settings": {
                    "solving_technique": "normal_equations"                              // 'normal_equations', 'qr_decomposition'
                }
            },
            "HROM":{
                "element_selection_type": "empirical_cubature",
                "element_selection_svd_truncation_tolerance": 1e-30,
                "create_hrom_visualization_model_part" : true,
                "echo_level" : 0
            }
        }""")
    return general_rom_manager_parameters

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# mult params
#
def get_multiple_params_by_Halton_test(number_of_values,angle,mach):
    sampler = qmc.Halton(d=2)
    sample = sampler.random(number_of_values)
    mu = []
    values = qmc.scale(sample, [angle[0],mach[0]], [angle[1],mach[1]])
    for i in range(number_of_values):
        #Angle of attack , Mach infinit, Upwind factor constant
        mu.append([np.round(values[i,0],3), np.round(values[i,1],3), np.round(1.000,3), i, "test"])
        GenerateMeshe(i, np.round(values[i,0],3), "test")
    return mu

def get_multiple_params_by_Halton_train(number_of_values,angle,mach):
    sampler = qmc.Halton(d=2)
    sample = sampler.random(number_of_values)
    mu = []
    values = qmc.scale(sample, [angle[0],mach[0]], [angle[1],mach[1]])
    values[0,0] = angle[0]
    values[0,1] = mach[0]
    values[1,0] = angle[0]
    values[1,1] = mach[1]
    values[number_of_values-1,0] = angle[1]
    values[number_of_values-1,1] = mach[1]
    values[number_of_values-2,0] = angle[1]
    values[number_of_values-2,1] = mach[0]
    for i in range(number_of_values):
        #Angle of attack , Mach infinit
        mu.append([np.round(values[i,0],3), np.round(values[i,1],3), np.round(1.000,3), i, "train"])
        GenerateMeshe(i, np.round(values[i,0],3), "train")
    return mu

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def plot_mu_values(mu_train,mu_test):
    mu_train_a = np.zeros(len(mu_train))
    mu_train_m = np.zeros(len(mu_train))
    mu_test_a  = np.zeros(len(mu_test))
    mu_test_m  = np.zeros(len(mu_test))
    for i in range(len(mu_train)):
        mu_train_a[i] = mu_train[i][0]
        mu_train_m[i] = mu_train[i][1]
    for i in range(len(mu_test)):
        mu_test_a[i] = mu_test[i][0]
        mu_test_m[i] = mu_test[i][1]
    plt.plot(mu_train_m, mu_train_a, 'bs', label="Train Values")
    plt.plot(mu_test_m, mu_test_a, 'ro', label="Test Values")
    plt.title('Mu Values')
    plt.ylabel('Alpha')
    plt.xlabel('Mach')
    plt.grid(True)
    plt.show(block=False)
    plt.legend(bbox_to_anchor=(.85, 1.03, 1., .102), loc='upper left', borderaxespad=0.)
    plt.savefig("MuValues.png")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

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
    
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# plot
# 
def plot_Cps(mu_train,mu_test):
    case_names = ["FOM","ROM","HROM"]
    markercolor = ["ob","xr","+g"]
    for j in range(len(mu_train)):
        cp_min = 0
        cp_max = 0
        fig = plt.figure()
        fig.set_figwidth(12.0)
        fig.set_figheight(8.0)
        for n, name in enumerate(case_names):
            if os.path.exists("Data/"+name+"_"+"train"+"_"+str(j)+".dat"):
                x  = np.loadtxt("Data/"+name+"_"+"train"+"_"+str(j)+".dat",usecols=(0,))
                cp = np.loadtxt("Data/"+name+"_"+"train"+"_"+str(j)+".dat",usecols=(3,))
                fig = plt.plot(x, cp, markercolor[n], markersize = 3.0, label = name)
                if np.min(cp) < cp_min:
                    cp_min = np.min(cp)
                if np.max(cp) > cp_max:
                    cp_max = np.max(cp)
        fig = plt.title('Cp vs x - ' + "angle: " + str(mu_train[j][0]) + "º " + "mach: " + str(mu_train[j][1]))
        fig = plt.axis([-0.05,1.35,cp_max+0.1,cp_min-0.1])
        fig = plt.ylabel('Cp')
        fig = plt.xlabel('x')
        fig = plt.grid()
        fig = plt.legend()
        fig = plt.tight_layout()
        fig = plt.savefig("Captures/Simulation_train_" + str(j) + "_A_" + str(mu_train[j][0]) + "_M_" + str(mu_train[j][1])+".png")
        fig = plt.close('all')

    for j in range(len(mu_test)):
        cp_min = 0
        cp_max = 0
        fig = plt.figure()
        fig.set_figwidth(12.0)
        fig.set_figheight(8.0)
        for n, name in enumerate(case_names):
            if os.path.exists("Data/"+name+"_"+"test"+"_"+str(j)+".dat"):
                x  = np.loadtxt("Data/"+name+"_"+"test"+"_"+str(j)+".dat",usecols=(0,))
                cp = np.loadtxt("Data/"+name+"_"+"test"+"_"+str(j)+".dat",usecols=(3,))
                fig = plt.plot(x, cp, markercolor[n], markersize = 3.0, label = name)
                if np.min(cp) < cp_min:
                    cp_min = np.min(cp)
                if np.max(cp) > cp_max:
                    cp_max = np.max(cp)
        fig = plt.title('Cp vs x - ' + "angle: " + str(mu_test[j][0]) + "º " + "mach: " + str(mu_test[j][1]))
        fig = plt.axis([-0.05,1.35,cp_max+0.1,cp_min-0.1])
        fig = plt.ylabel('Cp')
        fig = plt.xlabel('x')
        fig = plt.grid()
        fig = plt.legend()
        fig = plt.tight_layout()
        fig = plt.savefig("Captures/Simulation_test_" + str(j) + "_A_" + str(mu_test[j][0]) + "_M_" + str(mu_test[j][1])+".png")
        fig = plt.close('all')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# save / load
# 
def save_mu_parameters(mu_train,mu_test):
    archivo = open('Data/mu_train.dat', 'wb')
    pickle.dump(mu_train, archivo)
    archivo.close()
    archivo = open('Data/mu_test.dat', 'wb')
    pickle.dump(mu_test, archivo)
    archivo.close()

def load_mu_parameters():
    archivo = open('Data/mu_train.dat', 'rb')
    mu_train = pickle.load(archivo)
    archivo.close()
    archivo = open('Data/mu_test.dat', 'rb')
    mu_test = pickle.load(archivo)
    archivo.close()
    return mu_train,mu_test

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
def CleanFolder():
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Results')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Data')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Captures')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Meshes')
    # KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('')
    os.mkdir("Results")
    os.mkdir("Data")
    os.mkdir("Captures")
    os.mkdir("Meshes")
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


if __name__ == "__main__":

    # Minimo 5 por las cuatro esquinas y un punto interno para entrenar
    NumberofMuTrain = 40
    NumberOfMuTest  = 40

    load_old_mu_parameters = False

    # Definir rango de valores de mach y angulo de ataque
    mach_range  = [0.60, 0.75]
    angle_range = [0.00, 2.00] 

    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    if load_old_mu_parameters:
        mu_train, mu_test = load_mu_parameters()
    else:
        CleanFolder()
        
        mu_train = get_multiple_params_by_Halton_train(NumberofMuTrain,angle_range,mach_range) 
        mu_test  = get_multiple_params_by_Halton_test(NumberOfMuTest,angle_range,mach_range) 

        save_mu_parameters(mu_train,mu_test)
        plot_mu_values(mu_train,mu_test)

    general_rom_manager_parameters = GetRomManagerParameters()

    project_parameters_name = "ProjectParametersPrimalROM.json"

    rom_manager = RomManager(project_parameters_name,general_rom_manager_parameters,CustomizeSimulation,UpdateProjectParameters, UpdateMaterialParametersFile)

    rom_manager.Fit(mu_train)

    rom_manager.Test(mu_test)

    rom_manager.PrintErrors()

    plot_Cps(mu_train,mu_test)


