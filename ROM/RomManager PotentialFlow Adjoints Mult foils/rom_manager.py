import KratosMultiphysics
import numpy as np
from scipy.stats import qmc
import matplotlib.pyplot as plt
import math
import KratosMultiphysics.kratos_utilities
from KratosMultiphysics.MeshMovingApplication.mesh_moving_analysis import MeshMovingAnalysis
from KratosMultiphysics.RomApplication.rom_manager import RomManager

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def CustomizeSimulation(cls, global_model, parameters):

    class CustomSimulation(cls):

        def __init__(self, model,project_parameters, custom_param = None):
            super().__init__(model,project_parameters)
            self.custom_param  = custom_param

        def Initialize(self):

            angle_of_attack = parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].GetDouble()
            parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(0.0) 
            with open("ProjectParametersMeshMoving.json",'r') as parameter_file:
                mesh_parameters = KratosMultiphysics.Parameters(parameter_file.read()) 
            mesh_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["rotation_angle"].SetDouble(angle_of_attack)
            mesh_simulation = MeshMovingAnalysis(self.model,mesh_parameters)

            super().Initialize()

            mesh_simulation.Run()

        def FinalizeSolutionStep(self):
            super().FinalizeSolutionStep()

        def CustomMethod(self):
            return self.custom_param

    return CustomSimulation(global_model, parameters)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def UpdateProjectParameters(parameters, mu=None, Id=None, name=None):
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(mu[0])
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mu[1])
    parameters["processes"]["list_other_processes"][0]["Parameters"]["file_settings"]["file_name"].SetString(name+"_Primal_Data/MainModelPart_"+str(Id))
    return parameters

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def GetRomManagerParameters():
    general_rom_manager_parameters = KratosMultiphysics.Parameters("""{
            "rom_stages_to_train" : ["ROM","HROM"],      // ["ROM","HROM"]
            "rom_stages_to_test"  : ["ROM","HROM"],      // ["ROM","HROM"]
            "paralellism" : null,                        // null, TODO: add "compss"
            "projection_strategy": "galerkin",           // "lspg", "galerkin", "petrov_galerkin"
            "save_gid_output": false,                     // false, true #if true, it must exits previously in the ProjectParameters.json
            "save_vtk_output": false,                    // false, true #if true, it must exits previously in the ProjectParameters.json
            "output_name": "id",                         // "id" , "mu"
            "ROM":{
                "svd_truncation_tolerance": 1e-6,
                "model_part_name": "MainModelPart",                                      // This changes depending on the simulation: Structure, FluidModelPart, ThermalPart #TODO: Idenfity it automatically
                "nodal_unknowns": [ "VELOCITY_POTENTIAL",
                                    "AUXILIARY_VELOCITY_POTENTIAL"],                     // Main unknowns. Snapshots are taken from these
                "rom_basis_output_format": "json",                                       // "json" "numpy"
                "rom_basis_output_name": "PrimalRomParameters",
                "snapshots_control_type": "step",                                        // "step", "time"
                "snapshots_interval": 1,
                "petrov_galerkin_training_parameters":{
                    "basis_strategy": "residuals",                                        // 'residuals', 'jacobian'
                    "include_phi": false,
                    "svd_truncation_tolerance": 1e-9,
                    "echo_level": 0
                }
            },
            "HROM":{
                "element_selection_type": "empirical_cubature",
                "element_selection_svd_truncation_tolerance": 1e-9,
                "create_hrom_visualization_model_part" : false,
                "echo_level" : 0
            }
        }""")
    return general_rom_manager_parameters

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# mult params
#    
def get_multiple_params_by_Halton_test(number_of_values):
    sampler = qmc.Halton(d=2)
    sample = sampler.random(number_of_values)
    #Angle of attack
    l_angle = -1.0
    u_angle =  6.0
    #Mach infinit
    l_mach = 0.03
    u_mach = 0.6
    mu = []
    values = qmc.scale(sample, [l_angle,l_mach], [u_angle,u_mach])
    for i in range(number_of_values):
        #Angle of attack , Mach infinit
        mu.append([values[i,0] * math.pi / 180.0, values[i,1]])
    return mu

def get_multiple_params_by_Halton_train(number_of_values):
    sampler = qmc.Halton(d=2)
    sample = sampler.random(number_of_values)
    #Angle of attack
    l_angle = -1.0
    u_angle =  6.0
    #Mach infinit
    l_mach = 0.03
    u_mach = 0.6
    mu = []
    values = qmc.scale(sample, [l_angle,l_mach], [u_angle,u_mach])
    values[0,0] = l_angle
    values[0,1] = l_mach
    values[1,0] = l_angle
    values[1,1] = u_mach
    values[number_of_values-1,0] = u_angle
    values[number_of_values-1,1] = u_mach
    values[number_of_values-2,0] = u_angle
    values[number_of_values-2,1] = l_mach
    for i in range(number_of_values):
        #Angle of attack , Mach infinit
        mu.append([values[i,0] * math.pi / 180.0, values[i,1]])
    return mu
    
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def save_mu_parameters(mu_train,mu_test):
    import pickle
    archivo = open('mu_train.dat', 'wb')
    pickle.dump(mu_train, archivo)
    archivo.close()
    archivo = open('mu_test.dat', 'wb')
    pickle.dump(mu_test, archivo)
    archivo.close()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def load_mu_parameters():
    import pickle
    archivo = open('mu_train.dat', 'rb')
    lista = pickle.load(archivo)
    mu_train = np.asarray(lista)
    archivo.close()
    archivo = open('mu_test.dat', 'rb')
    lista = pickle.load(archivo)
    mu_test = np.asarray(lista)
    archivo.close()
    return mu_train,mu_test
    

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def plot_mu_values(mu_train,mu_test,name):
    mu_train_a = np.zeros(len(mu_train))
    mu_train_m = np.zeros(len(mu_train))
    mu_test_a  = np.zeros(len(mu_test))
    mu_test_m  = np.zeros(len(mu_test))
    for i in range(len(mu_train)):
        mu_train_a[i] = -mu_train[i][0] * 180 / math.pi + 5.0
        mu_train_m[i] = mu_train[i][1]
    for i in range(len(mu_test)):
        mu_test_a[i] = -mu_test[i][0] * 180 / math.pi + 5.0
        mu_test_m[i] = mu_test[i][1]
    plt.plot(mu_train_m, mu_train_a, 'bs', label="Train Values")
    plt.plot(mu_test_m, mu_test_a, 'ro', label="Test Values")
    plt.title('Mu Values')
    plt.ylabel('Alpha')
    plt.xlabel('Mach')
    plt.grid()
    # plt.show(block=False)
    # plt.tight_layout()
    plt.legend()
    plt.legend(bbox_to_anchor=(.85, 1.03, 1., .102), loc='upper left', borderaxespad=0.)
    plt.savefig(name)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


if __name__ == "__main__":
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Results')

    mu_train = get_multiple_params_by_Halton_train(25) 
    mu_test  = get_multiple_params_by_Halton_test(25) 

    save_mu_parameters(mu_train,mu_test)

    # mu_train, mu_test = load_mu_parameters()

    plot_mu_values(mu_train,mu_test,"Mu_Values.png")

    general_rom_manager_parameters = GetRomManagerParameters()

    primal_project_parameters_name = "ProjectParametersPrimalROM.json"

    primal_rom_manager = RomManager(primal_project_parameters_name,general_rom_manager_parameters,CustomizeSimulation,UpdateProjectParameters,"Primal_")

    primal_rom_manager.Fit(mu_train) 

    primal_rom_manager.Test(mu_test) 

    primal_rom_manager.PrintErrors()


    adjoint_project_parameters_name = "ProjectParametersAdjointROM.json"

    general_rom_manager_parameters["ROM"]["nodal_unknowns"].SetStringArray([ "ADJOINT_VELOCITY_POTENTIAL"
                                                                            ,"ADJOINT_AUXILIARY_VELOCITY_POTENTIAL"])

    general_rom_manager_parameters["ROM"]["rom_basis_output_name"].SetString("AdjointRomParameters")

    adjoint_rom_manager = RomManager(adjoint_project_parameters_name,general_rom_manager_parameters,CustomizeSimulation,UpdateProjectParameters,"Adjoint_")

    adjoint_rom_manager.Fit(mu_train) 

    adjoint_rom_manager.Test(mu_test) 

    adjoint_rom_manager.PrintErrors()

