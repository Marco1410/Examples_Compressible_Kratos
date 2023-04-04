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
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(mu[0])
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mu[1])

    return parameters

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def GetRomManagerParameters():
    """
    This function allows to easily modify all the parameters for the ROM simulation.
    The returned KratosParameter object is seamlessly used inside the RomManager.
    """
    general_rom_manager_parameters = KratosMultiphysics.Parameters("""{
            "rom_stages_to_train" : ["ROM"],      // ["ROM","HROM"]
            "rom_stages_to_test"  : ["ROM"],      // ["ROM","HROM"]
            "paralellism" : null,                        // null, TODO: add "compss"
            "projection_strategy": "lspg",           // "lspg", "galerkin", "petrov_galerkin"
            "save_gid_output": true,                     // false, true #if true, it must exits previously in the ProjectParameters.json
            "save_vtk_output": false,                    // false, true #if true, it must exits previously in the ProjectParameters.json
            "output_name": "mu",                         // "id" , "mu"
            "ROM":{
                "svd_truncation_tolerance": 1e-6,
                "model_part_name": "MainModelPart",                                      // This changes depending on the simulation: Structure, FluidModelPart, ThermalPart #TODO: Idenfity it automatically
                "nodal_unknowns": ["VELOCITY_POTENTIAL","AUXILIARY_VELOCITY_POTENTIAL"], // Main unknowns. Snapshots are taken from these
                "rom_basis_output_format": "json",                                      //TODO: add "numpy"
                "rom_basis_output_name": "RomParameters",
                "snapshots_control_type": "step",                                        // "step", "time"
                "snapshots_interval": 1,
                "petrov_galerkin_training_parameters":{
                    "basis_strategy": "residuals",                                        // 'residuals', 'jacobian'
                    "include_phi": false,
                    "svd_truncation_tolerance": 1e-6,
                    "echo_level": 0
                }
            },
            "HROM":{
                "element_selection_type": "empirical_cubature",
                "element_selection_svd_truncation_tolerance": 1e-6,
                "create_hrom_visualization_model_part" : false,
                "echo_level" : 0
            }
        }""")

    return general_rom_manager_parameters

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

if __name__ == "__main__":
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Results')

    def random_samples_from_interval(initial, final, number_of_samples):
        import numpy as np
        return initial + np.random.rand(number_of_samples)*(final-initial)

    def get_multiple_params_test():
        number_of_params_1 = 1
        number_of_params_2 = 3
        param1 = random_samples_from_interval(-3.0, -3.0,number_of_params_1)
        param2 = random_samples_from_interval( 0.03, 0.6,number_of_params_2)
        mu = []
        for i in range(number_of_params_1):
            for j in range(number_of_params_2):
                mu.append([param1[i] * math.pi / 180.0, param2[j]])
        return mu

    def get_multiple_params_train():
        number_of_values = 4
        sampler = qmc.Halton(d=1)
        sampler.seed = 10
        sample = sampler.random(number_of_values)
        #Mach infinit
        l_mach = [0.03]
        u_mach = [0.6]
        mu = []
        values = qmc.scale(sample, [l_mach[0]], [u_mach[0]])
        values[0] = l_mach
        values[number_of_values-1] = u_mach
        for i in range(number_of_values):
            mu.append([-3.0 * math.pi / 180.0, values[i]])
        return mu
    
    # def get_multiple_params_test():
    #     number_of_params_1 = 1
    #     number_of_params_2 = 1
    #     param1 = random_samples_from_interval(-6.0, 1.0,number_of_params_1)
    #     param2 = random_samples_from_interval( 0.3, 0.3,number_of_params_2)
    #     mu = []
    #     for i in range(number_of_params_1):
    #         for j in range(number_of_params_2):
    #             mu.append([param1[i] * math.pi / 180.0, param2[j]])
    #     return mu
    
    # def get_multiple_params_train():
    #     #Angle of attack
    #     mu = []
    #     #mu.append([-6.0 * math.pi / 180.0, 0.3])
    #     #mu.append([-5.5 * math.pi / 180.0, 0.3])
    #    # mu.append([-5.0 * math.pi / 180.0, 0.3])
    #    # mu.append([-4.5 * math.pi / 180.0, 0.3])
    #     #mu.append([-4.0 * math.pi / 180.0, 0.3])
    #     #mu.append([-3.5 * math.pi / 180.0, 0.3])
    #     #mu.append([-3.0 * math.pi / 180.0, 0.3])
    #    # mu.append([-2.5 * math.pi / 180.0, 0.3])
    #     mu.append([-2.0 * math.pi / 180.0, 0.3])
    #    # mu.append([-1.5 * math.pi / 180.0, 0.3])
    #     #mu.append([-1.0 * math.pi / 180.0, 0.3])
    #     #mu.append([-0.5 * math.pi / 180.0, 0.3])
    #     #mu.append([ 0.0 * math.pi / 180.0, 0.3])
    #     #mu.append([ 0.5 * math.pi / 180.0, 0.3])
    #     #mu.append([ 1.0 * math.pi / 180.0, 0.3])
    #     return mu

    mu_train = get_multiple_params_train() # random train parameters

    mu_test  = get_multiple_params_test() #random test parameters

    plot_mu_values = True

    if plot_mu_values:
        mu_train_a = np.zeros(len(mu_train))
        mu_train_m = np.zeros(len(mu_train))
        mu_test_a  = np.zeros(len(mu_test))
        mu_test_m  = np.zeros(len(mu_test))
        for i in range(len(mu_train)):
            mu_train_a[i] = mu_train[i][0] * 180/math.pi + 5.0
            mu_train_m[i] = mu_train[i][1]
        for i in range(len(mu_test)):
            mu_test_a[i] = mu_test[i][0] * 180/math.pi + 5.0
            mu_test_m[i] = mu_test[i][1]
        plt.plot(mu_train_m, mu_train_a, 'bs', label="Train Values")
        plt.plot(mu_test_m, mu_test_a, 'ro', label="Test Values")
        plt.title('Mu Values')
        plt.ylabel('Alpha')
        plt.xlabel('Mach')
        plt.grid(True)
        plt.show(block=False)
        plt.legend(loc='best')
        plt.savefig("MuValues.png")


    general_rom_manager_parameters = GetRomManagerParameters()

    project_parameters_name = "ProjectParametersPrimalROM.json"

    rom_manager = RomManager(project_parameters_name,general_rom_manager_parameters,CustomizeSimulation,UpdateProjectParameters)

    """if no list "mu" is passed, the case already contained in the ProjectParametes and CustomSimulation is launched (useful for example for a single time dependent simulation)"""

    rom_manager.Fit(mu_train) 

    rom_manager.Test(mu_test) 

    rom_manager.PrintErrors()


