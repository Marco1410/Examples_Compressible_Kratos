import os
import time
import openpyxl
import KratosMultiphysics.kratos_utilities
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
import KratosMultiphysics
from local_rom_manager import RomManager 
from parameters import * 
from Fom_database import * 
from PEBL_RBF import * 
from KratosMultiphysics.RomApplication.calculate_rom_basis_output_process import CalculateRomBasisOutputProcess


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


def CustomizeSimulation(cls, global_model, parameters):

    class CustomSimulation(cls):

        def __init__(self, model,project_parameters):
            super().__init__(model,project_parameters)

        def Run(self):
            if self._GetSimulationName() == "Analysis": # FOM
                self.Initialize()
            elif self._GetSimulationName() == "::[ROM Simulation]:: ": # ROM

                start_time = time.time()
                self.Initialize()
                self.RunSolutionLoop()
                self.Finalize()
                exe_time = time.time() - start_time

                if parameters["output_processes"].Has("gid_output"):

                    simulation_name = parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].GetString().removeprefix('Results/')
                    cluster_path    = parameters["output_processes"]["rom_output"][0]['Parameters']['rom_basis_output_folder'].GetString()
                    cluster_number  = cluster_path.removeprefix('rom_bases/rom_data_cluster_')

                    angle = parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].GetDouble()
                    mach  = parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].GetDouble()

                    fom_name = f'{angle}, {mach}'

                    info_steps_list = []
                    if 'HROM' in simulation_name:
                        if 'Fit' in simulation_name:
                            data_output_name = 'hrom_data.xlsx'
                        elif 'Test' in simulation_name:
                            data_output_name = 'test_hrom_data.xlsx'
                    elif 'ROM' in simulation_name:
                        if 'Fit' in simulation_name:
                            data_output_name = 'rom_data.xlsx'
                        elif 'Test' in simulation_name:
                            data_output_name = 'test_rom_data.xlsx' 

                    for process in self._GetListOfOutputProcesses():
                            if isinstance(process, CalculateRomBasisOutputProcess):
                                BasisOutputProcess = process
                    fom = np.load(f'DataBase/Snapshots/{fom_name}.npy')
                    rom = BasisOutputProcess._GetSnapshotsMatrix()

                    modes = np.load(f'{cluster_path}/RightBasisMatrix.npy').shape[1]
                
                    info_steps_list.append([cluster_number,
                                            angle,
                                            mach, 
                                            self.model["MainModelPart"].ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER],
                                            self.model["MainModelPart"].ProcessInfo[KratosMultiphysics.RESIDUAL_NORM],
                                            np.linalg.norm(fom-rom)/np.linalg.norm(fom),
                                            modes,
                                            round(exe_time, 2)])
                    
                    if os.path.exists(data_output_name):
                        wb = openpyxl.load_workbook(data_output_name)
                        hoja = wb.active
                        for item in info_steps_list:
                            hoja.append(item)
                        wb.save(data_output_name)
                    else:
                        wb = openpyxl.Workbook()
                        hoja = wb.active
                        hoja.append(('Cluster', 'Angle [º]', 'Mach', 'NL iterations', 'Residual norm', 'Approximation error [%]', 'Modes', 'Time [sec]'))
                        for item in info_steps_list:
                            hoja.append(item)
                        wb.save(data_output_name)

                    if not os.path.exists(f'Skin_Data'):
                        os.mkdir(f'Skin_Data')
                    if not os.path.exists(f'Captures'):
                        os.mkdir(f'Captures')
                    skin_data_filename = f"Skin_Data/{simulation_name}.dat"
                    capture_filename   = f"Captures/{simulation_name.split('_')[1]}.png"

                    fout = open(skin_data_filename,'w')
                    modelpart = self.model["MainModelPart.Body2D_Body"]
                    for node in modelpart.Nodes:
                        x = node.X ; y = node.Y ; z = node.Z
                        cp = node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
                        fout.write("%s %s %s %s\n" %(x,y,z,cp))
                    fout.close()

                    #### ONLINE CP PLOT
                    ######################################################################
                    cp_min = 0
                    cp_max = 0
                    fig = plt.figure()
                    fig.set_figwidth(12.0)
                    fig.set_figheight(8.0)
                    #### FOM ######
                    fom_skin_data_filename = f"DataBase/FOM_Skin_Data/{fom_name}.dat"
                    x_fom  = np.loadtxt(fom_skin_data_filename, usecols=(0,))
                    cp_fom = np.loadtxt(fom_skin_data_filename, usecols=(3,))
                    fig = plt.plot(x_fom, cp_fom, 'ob', markersize = 2.0, label = 'FOM')
                    #### ROM ######
                    if 'HROM' in simulation_name:
                        rom_skin_data_filename = f"Skin_Data/ROM_{simulation_name.split('_')[1]}.dat"
                        x_rom  = np.loadtxt(rom_skin_data_filename, usecols=(0,))
                        cp_rom = np.loadtxt(rom_skin_data_filename, usecols=(3,))
                        fig = plt.plot(x_rom, cp_rom, 'xr', markersize = 2.0, label = 'ROM')
                    ###############
                    x  = np.loadtxt(skin_data_filename, usecols=(0,))
                    cp = np.loadtxt(skin_data_filename, usecols=(3,))
                    fig = plt.plot(x, cp, '+g', markersize = 2.0, label = simulation_name.split('_')[0])
                    if np.min(cp) < cp_min:
                        cp_min = np.min(cp)
                    if np.max(cp) > cp_max:
                        cp_max = np.max(cp)
                    fig = plt.title('Cp vs x')
                    fig = plt.axis([-0.05,1.35,cp_max+0.1,cp_min-0.1])
                    fig = plt.ylabel('Cp')
                    fig = plt.xlabel('x')
                    fig = plt.grid()
                    fig = plt.legend()
                    fig = plt.tight_layout()
                    fig = plt.savefig(capture_filename)
                    fig = plt.close('all')

    return CustomSimulation(global_model, parameters)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def UpdateProjectParameters(parameters, mu=None):
    angle_of_attack        = mu[0]
    mach_infinity          = mu[1]
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(angle_of_attack)
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach_infinity)
    return parameters

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def UpdateMaterialParametersFile(material_parametrs_file_name, mu):
    pass

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


def GetRomManagerParameters():
    general_rom_manager_parameters = KratosMultiphysics.Parameters("""{
            "rom_stages_to_train" : ["ROM"],             // ["ROM","HROM"]
            "rom_stages_to_test"  : ["ROM"],             // ["ROM","HROM"]
            "paralellism" : null,                        // null, TODO: add "compss"
            "projection_strategy": "lspg",            // "lspg", "galerkin", "petrov_galerkin"
            "assembling_strategy": "global",            // "global", "elemental"
            "save_gid_output": true,                    // false, true #if true, it must exits previously in the ProjectParameters.json
            "save_vtk_output": false,                    // false, true #if true, it must exits previously in the ProjectParameters.json
            "output_name": "mu",                         // "id" , "mu"
            "ROM":{
                "svd_truncation_tolerance": 1e-12,
                "model_part_name": "MainModelPart",                            // This changes depending on the simulation: Structure, FluidModelPart, ThermalPart #TODO: Idenfity it automatically
                "nodal_unknowns": ["VELOCITY_POTENTIAL","AUXILIARY_VELOCITY_POTENTIAL"],     // Main unknowns. Snapshots are taken from these
                "rom_basis_output_format": "numpy",
                "rom_basis_output_name": "RomParameters",
                "rom_basis_output_folder": "rom_data",
                "snapshots_control_type": "step",                          // "step", "time"
                "snapshots_interval": 1,
                "galerkin_rom_bns_settings": {
                    "monotonicity_preserving": false
                },
                "lspg_rom_bns_settings": {
                    "train_petrov_galerkin": false,
                    "basis_strategy": "reactions",                        // 'residuals', 'jacobian', 'reactions'
                    "include_phi": false,
                    "svd_truncation_tolerance": 1e-12,
                    "solving_technique": "normal_equations",              // 'normal_equations', 'qr_decomposition'
                    "monotonicity_preserving": false
                },
                "petrov_galerkin_rom_bns_settings": {
                    "monotonicity_preserving": false
                }
            },
            "HROM":{
                "element_selection_type": "empirical_cubature",
                "element_selection_svd_truncation_tolerance": 1e-12,
                "create_hrom_visualization_model_part" : false,
                "echo_level" : 0
            }
        }""")

    return general_rom_manager_parameters


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


if __name__ == "__main__":

    ###############################
    # PARAMETERS SETTINGS
    number_of_mu_train = 100
    mach_range         = [ 0.70, 0.75]
    angle_range        = [ 1.00, 2.00]
    ###############################

    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    sampler = qmc.Halton(d=2)
    mu_train = []

    if number_of_mu_train > 0:
        sample = sampler.random(number_of_mu_train)
        values = qmc.scale(sample, [angle_range[0],mach_range[0]], [angle_range[1],mach_range[1]])
        for i in range(number_of_mu_train):
            #Angle of attack , Mach infinit
            mu_train.append([values[i,0], values[i,1]])

        np.save(f'mu_train', mu_train)

    plot_mu_values(mu_train, [], 'MuValues')

    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    os.mkdir(f"DataBase")
    os.mkdir(f"DataBase/Snapshots")
    os.mkdir(f"DataBase/Snapshots_not_converged")
    os.mkdir(f"DataBase/FOM_Results")
    os.mkdir(f"DataBase/FOM_Skin_Data")
    os.mkdir(f"DataBase/FOM_Captures")

    LaunchFOMSimulations(list(mu_train)) 

    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    general_rom_manager_parameters = GetRomManagerParameters()
    project_parameters_name = "ProjectParameters.json"
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    rom_manager = RomManager(   project_parameters_name, 
                                general_rom_manager_parameters, 
                                CustomizeSimulation, 
                                UpdateProjectParameters, 
                                UpdateMaterialParametersFile )          
        
    old_fom_rom_train_error,_ = rom_manager.Fit(mu_train)
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    if old_fom_rom_train_error > 1e-12:
        iter = 0 
        os.mkdir(f"Mu_parameters")
        np.save(f'Mu_parameters/mu_{iter}', mu_train)
        plot_mu_values(mu_train, [], f'Mu_parameters/Mu_{iter}')

        new_mu = mu_train
        new_fom_rom_train_error = 1.0

        while new_fom_rom_train_error > 1e-9:
            iter += 1
            KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('rom_data.xlsx')  
            #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            new_value = []
            sample = sampler.random(5)
            values = qmc.scale(sample, [angle_range[0],mach_range[0]], [angle_range[1],mach_range[1]])
            for i in range(len(values)):
            #Angle of attack , Mach infinit
                new_value.append([values[i,0], values[i,1]])

            LaunchFOMSimulations(new_value)
            #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::     
            new_fom_rom_train_error,_ = rom_manager.Fit(new_mu + new_value)
            #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            if new_fom_rom_train_error < old_fom_rom_train_error:
                old_fom_rom_train_error = new_fom_rom_train_error
                new_mu = new_mu + new_value
                mu_train = new_mu
                np.save(f'Mu_parameters/mu_{iter}', new_mu)
                plot_mu_values(new_mu, [], f'Mu_parameters/Mu_{iter}')
        
    modes = np.load(f'rom_data/RightBasisMatrix.npy').shape[1]
    fom_rom_train_errors = []           
    fom_rom_train_errors.append([f'Train error FOM - ROM, Cluster {i}: {old_fom_rom_train_error:.2E}, Modes: {modes}, Points: {len(mu_train)}'])


    print('::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
    for i in range(len(fom_rom_train_errors)):
        print(fom_rom_train_errors[i])  
    print('::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
    
    np.save('fom_rom_train_errors', fom_rom_train_errors)
