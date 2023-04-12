import KratosMultiphysics
import numpy as np
from scipy.stats import qmc
import matplotlib.pyplot as plt
import math
import KratosMultiphysics.kratos_utilities
from KratosMultiphysics.RomApplication.rom_manager import RomManager

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def CustomizeSimulation(cls, global_model, parameters):

    class CustomSimulation(cls):

        def __init__(self, model,project_parameters, custom_param = None):
            super().__init__(model,project_parameters)
            self.custom_param  = custom_param
            """
            Customize as needed
            """

        def Initialize(self):
            super().Initialize()
            """
            Customize as needed
            """

        def FinalizeSolutionStep(self):
            super().FinalizeSolutionStep()
            """
            Customize as needed
            """


        def CustomMethod(self):
            """
            Customize as needed
            """
            return self.custom_param

    return CustomSimulation(global_model, parameters)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def UpdateProjectParameters(parameters, mu=None):
    """
    Customize ProjectParameters here for imposing different conditions to the simulations as needed
    """
    parameters["processes"]["ale_bc_process_list"][0]["Parameters"]["rotation_angle"].SetDouble(mu[0])
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mu[1])

    return parameters

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def GetRomManagerParameters():
    """
    This function allows to easily modify all the parameters for the ROM simulation.
    The returned KratosParameter object is seamlessly used inside the RomManager.
    """
    general_rom_manager_parameters = KratosMultiphysics.Parameters("""{
            "rom_stages_to_train" : [""],      // ["ROM","HROM"]
            "rom_stages_to_test"  : ["ROM"],      // ["ROM","HROM"]
            "projection_strategy": "galerkin",           // "lspg", "galerkin", "petrov_galerkin"
            "save_gid_output": false,                     // false, true #if true, it must exits previously in the ProjectParameters.json
            "save_vtk_output": false,                    // false, true #if true, it must exits previously in the ProjectParameters.json
            "output_name": "id",                         // "id" , "mu"
            "ROM":{
                "svd_truncation_tolerance": 1e-3,
                "model_part_name": "MainModelPart",                                      // This changes depending on the simulation: Structure, FluidModelPart, ThermalPart #TODO: Idenfity it automatically
                "nodal_unknowns": ["VELOCITY_POTENTIAL","AUXILIARY_VELOCITY_POTENTIAL"], // Main unknowns. Snapshots are taken from these
                "rom_basis_output_format": "json",                                      // "json" "numpy"
                "rom_basis_output_name": "RomParameters",
                "snapshots_control_type": "step",                                        // "step", "time"
                "snapshots_interval": 1,
                "petrov_galerkin_training_parameters":{
                    "basis_strategy": "residuals",                                        // 'residuals', 'jacobian'
                    "include_phi": false,
                    "svd_truncation_tolerance": 1e-3,
                    "echo_level": 0
                }
            },
            "HROM":{
                "element_selection_type": "empirical_cubature",
                "element_selection_svd_truncation_tolerance": 1e-3,
                "create_hrom_visualization_model_part" : false,
                "echo_level" : 0
            }
        }""")

    return general_rom_manager_parameters

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# mult params
#    
def get_multiple_params_random_test(number_of_params_1,number_of_params_2):
    #Angle of attack
    param1 = random_samples_from_interval(-1.0, 6.0,number_of_params_1)
    #Mach infinit
    param2 = random_samples_from_interval( 0.03, 0.6,number_of_params_2)
    mu = []
    for i in range(number_of_params_1):
        for j in range(number_of_params_2):
            #Angle of attack , Mach infinit
            mu.append([param1[i] * math.pi / 180.0, param2[j]])
    return mu

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

def get_multiple_params_by_angle_train(number_of_values):
    sampler = qmc.Halton(d=1)
    sample = sampler.random(number_of_values)
    #Angle of attack
    l_angle = -1.0
    u_angle =  6.0
    #Mach infinit
    l_mach = 0.03
    u_mach = 0.6
    mu = []

    params = np.zeros((15))
    for i in range(len(params)):
        params[i] = l_angle + i * 0.5
        
    for i in range(len(params)): 
        sample = sampler.random(number_of_values)
        values = qmc.scale(sample, [l_mach], [u_mach])
        #Angle of attack , Mach infinit
        mu.append([params[i] * math.pi / 180.0, l_mach]) 
        mu.append([params[i] * math.pi / 180.0, u_mach]) 
        for j in range(number_of_values):
            #Angle of attack , Mach infinit
            mu.append([params[i] * math.pi / 180.0, values[j]])  
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

# ## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# mult mach
#
def get_params_by_mach_test(number_of_params_2):
    #Angle of attack
    param1 = [3.0]
    #Mach infinit
    param2 = random_samples_from_interval( 0.03, 0.6,number_of_params_2)
    mu = []
    for i in range(number_of_params_2):
        #Angle of attack , Mach infinit
        mu.append([param1[0] * math.pi / 180.0, param2[i]])
    return mu

def get_params_by_mach_train(number_of_values):
    sampler = qmc.Halton(d=1)
    sample = sampler.random(number_of_values)
    #Mach infinit
    l_mach = 0.03
    u_mach = 0.6
    mu = []
    values = qmc.scale(sample, [l_mach], [u_mach])
    values[0] = l_mach
    values[number_of_values-1] = u_mach
    for i in range(number_of_values):
        #Angle of attack , Mach infinit
        mu.append([ 3.0 * math.pi / 180.0, values[i]])
    return mu

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# mult angle
#
def get_params_by_angle_test(number_of_params_1):
    #Angle of attack
    param1 = random_samples_from_interval(6.0, -1.0,number_of_params_1)
    #Mach infinit
    param2 = [0.3]
    mu = []
    for i in range(number_of_params_1):
        #Angle of attack , Mach infinit
        mu.append([param1[i] * math.pi / 180.0, param2[0]])
    return mu

def get_params_by_angle_train(number_of_values):
    sampler = qmc.Halton(d=1)
    sample = sampler.random(number_of_values)
    #Angle of attack
    l_angle =  6.0
    u_angle = -1.0
    mu = []
    values = qmc.scale(sample, [l_angle], [u_angle])
    values[0] = l_angle
    values[number_of_values-1] = u_angle
    for i in range(number_of_values):
        #Angle of attack , Mach infinit
        mu.append([values[i] * math.pi / 180.0, 0.3])
    return mu

def get_params_by_angle_fix_train():
    mu = []
    #Angle of attack , Mach infinit
    mu.append([6.0 * math.pi / 180.0, 0.3])
    #mu.append([5.75 * math.pi / 180.0, 0.3])
    #mu.append([5.5 * math.pi / 180.0, 0.3])
    #mu.append([5.25 * math.pi / 180.0, 0.3])
    mu.append([5.0 * math.pi / 180.0, 0.3])
    #mu.append([4.75 * math.pi / 180.0, 0.3])
    #mu.append([4.5 * math.pi / 180.0, 0.3])
    #mu.append([4.25 * math.pi / 180.0, 0.3])
    mu.append([4.0 * math.pi / 180.0, 0.3])
    #mu.append([3.75 * math.pi / 180.0, 0.3])
    #mu.append([3.5 * math.pi / 180.0, 0.3])
    #mu.append([3.25 * math.pi / 180.0, 0.3])
    mu.append([3.0 * math.pi / 180.0, 0.3])
    #mu.append([2.75 * math.pi / 180.0, 0.3])
    #mu.append([2.5 * math.pi / 180.0, 0.3])
    #mu.append([2.25 * math.pi / 180.0, 0.3])
    mu.append([2.0 * math.pi / 180.0, 0.3])
    #mu.append([1.75 * math.pi / 180.0, 0.3])
    #mu.append([1.5 * math.pi / 180.0, 0.3])
    #mu.append([1.25 * math.pi / 180.0, 0.3])
    mu.append([1.0 * math.pi / 180.0, 0.3])
    #mu.append([0.75 * math.pi / 180.0, 0.3])
    #mu.append([0.5 * math.pi / 180.0, 0.3])
    #mu.append([0.25 * math.pi / 180.0, 0.3])
    mu.append([ 0.0 * math.pi / 180.0, 0.3])
    #mu.append([-0.25 * math.pi / 180.0, 0.3])
    #mu.append([-0.5 * math.pi / 180.0, 0.3])
    #mu.append([-0.75 * math.pi / 180.0, 0.3])
    mu.append([-1.0 * math.pi / 180.0, 0.3])
    return mu

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def random_samples_from_interval(initial, final, number_of_samples):
    return initial + np.random.rand(number_of_samples)*(final-initial)
    
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def plot_mu_values(mu_train,mu_test):
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
    plt.grid(True)
    plt.show(block=False)
    plt.legend(bbox_to_anchor=(.85, 1.03, 1., .102), loc='upper left', borderaxespad=0.)
    plt.savefig("MuValues.png")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #




if __name__ == "__main__":
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Results')

    #################################################
    ############## GET PARAMETERS ###################
    #################################################

    ############## By angle and mach  ###############
    mu_train = get_multiple_params_by_Halton_train(5) 
    mu_test  = get_multiple_params_by_Halton_test(1) 

    # mu_train = get_multiple_params_by_Halton_train(5)
    # mu_test  = get_multiple_params_random_test(1,1)

    # mu_train = get_multiple_params_by_angle_train(2)
    # mu_test  = get_multiple_params_by_Halton_test(2) 

    ##############      By mach       ###############
    # mu_train = get_params_by_mach_train(3) 
    # mu_test  = get_params_by_mach_test(2)

    ##############      By angle      ###############
    # mu_train = get_params_by_angle_train(10)
    # mu_test  = get_params_by_angle_test(1)

    # mu_train = get_params_by_angle_fix_train()
    # mu_test  = get_params_by_angle_test(1) 
    #################################################

    plot_mu_values(mu_train,mu_test)

    general_rom_manager_parameters = GetRomManagerParameters()

    project_parameters_name = "ProjectParametersPrimalROM.json"

    rom_manager = RomManager(project_parameters_name,general_rom_manager_parameters,CustomizeSimulation,UpdateProjectParameters)

    rom_manager.Fit(mu_train) 

    rom_manager.Test(mu_test) 

    rom_manager.PrintErrors()


