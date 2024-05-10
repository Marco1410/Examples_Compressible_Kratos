import os
import time
import random
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
                    np.save(f'RomDataBase/{fom_name}',rom)

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
            "rom_stages_to_train" : ["ROM","HROM"],             // ["ROM","HROM"]
            "rom_stages_to_test"  : ["ROM","HROM"],             // ["ROM","HROM"]
            "paralellism" : null,                        // null, TODO: add "compss"
            "projection_strategy": "galerkin",            // "lspg", "galerkin", "petrov_galerkin"
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

    Update_parameters      = False
    Launch_fom_simulations = False
    Update_clusters        = True
    Train_rom              = True
    Test_rom               = True

    ###############################
    # PARAMETERS SETTINGS
    number_of_mu_train = 600
    number_of_mu_test  = 500
    mach_range         = [ 0.70, 0.75]
    angle_range        = [ 1.00, 2.00]
    ###############################
    # CLUSTERING SETTINGS
    bisection_tolerance = 0.95
    POD_tolerance       = 1e-12
    ###############################
    
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    if Update_parameters:
        KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_train.npy')
        KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_test.npy')
        KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_train_not_scaled.npy')
        KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_test_not_scaled.npy')
        KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('MuValues.png')
        KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('MuValuesNotScaled.png')

        mu_train, mu_test, mu_train_not_scaled, mu_test_not_scaled = get_multiple_parameters( number_train_values = number_of_mu_train,
                                                                                              number_test_values  = number_of_mu_test , 
                                                                                              angle               = angle_range       , 
                                                                                              mach                = mach_range        , 
                                                                                              method              = 'Halton'           )
    else:
        mu_train, mu_test, mu_train_not_scaled, mu_test_not_scaled = load_mu_parameters() 

    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    if Launch_fom_simulations:
        KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('DataBase')
        KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('fom_data.xlsx')
        os.mkdir(f"DataBase")
        os.mkdir(f"DataBase/Snapshots")
        os.mkdir(f"DataBase/Snapshots_not_converged")
        os.mkdir(f"DataBase/FOM_Results")
        os.mkdir(f"DataBase/FOM_Skin_Data")
        os.mkdir(f"DataBase/FOM_Captures")

        start_time = time.time()
        LaunchFOMSimulations(list(mu_train) + list(mu_test)) 
        exe_time = time.time() - start_time

        print(f' Executing took {round(exe_time, 2)} sec')

    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    if Update_clusters:
        KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Clustering_mu_test_plots')
        KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Mu_by_clusters')
        KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('predicted_indexes_list.npy')
        KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('PEBL clustering.png')
        os.mkdir(f"Clustering_mu_test_plots")
        os.mkdir(f"Mu_by_clusters")

        PEBL_Clustering( bisection_tolerance      = bisection_tolerance, 
                         POD_tolerance            = POD_tolerance, 
                         plot_custering_info      = True,
                         plot_clusters            = True,
                         mu_train                 = mu_train,
                         mu_train_not_scaled      = mu_train_not_scaled)

    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    fom_rom_train_errors = []; rom_hrom_train_errors = [];  fom_rom_test_errors = []; rom_hrom_test_errors = []

    # cluster_number = 0
    # Mu_by_clusters_list = ['mu_train_3']

    if Train_rom:
        KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('RomDataBase') 
        KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('RBF_interpolated_data')
        KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('rom_bases')
        KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Captures')
        KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Skin_Data')
        KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('fom_rom_train_errors.npy')
        KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('rom_hrom_train_errors.npy')
        KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('rom_data.xlsx')  
        KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('hrom_data.xlsx') 
        os.mkdir(f"RomDataBase")
        os.mkdir(f"RBF_interpolated_data")

        with os.scandir(f'Mu_by_clusters') as files:
            Mu_by_clusters_list = [file.name.removesuffix(".npy") for file in files if file.is_file() and file.name.endswith('.npy') and 'train' in file.name]
            Mu_by_clusters_list.sort()

        for i, mu_name in enumerate(Mu_by_clusters_list):
            mu = np.load(f'Mu_by_clusters/{mu_name}.npy')

            number_of_mu_by_cluster = 150
            if number_of_mu_by_cluster < len(mu):
                mu = np.array(random.sample(list(mu), number_of_mu_by_cluster))

            # RBF PREDICTION
            #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            mu_not_scaled = [mu_value for mu_value, mu_train_value in zip(mu_train_not_scaled, mu_train) if mu_train_value in mu]
            
            RBF_prediction(mu_train            = mu, 
                           mu_train_not_scaled = mu_not_scaled, 
                           mu_test             = mu,
                           mu_test_not_scaled  = mu_not_scaled)

            #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            general_rom_manager_parameters = GetRomManagerParameters()
            general_rom_manager_parameters['ROM']['rom_basis_output_folder'].SetString(f'rom_bases/rom_data_cluster_{i}')

            project_parameters_name = "ProjectParameters.json"

            rom_manager = RomManager(   project_parameters_name, 
                                        general_rom_manager_parameters, 
                                        CustomizeSimulation, 
                                        UpdateProjectParameters, 
                                        UpdateMaterialParametersFile )            
            #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

            fom_rom_train_error, rom_hrom_train_error = rom_manager.Fit(mu, store_all_snapshots = True)

            modes = np.load(f'rom_bases/rom_data_cluster_{i}/RightBasisMatrix.npy').shape[1]
            fom_rom_train_errors.append([f'Train error FOM - ROM, Cluster {i}: {fom_rom_train_error:.2E}, Modes: {modes}, Points: {len(mu)}'])
            if rom_hrom_train_error != 0:
                rom_hrom_train_errors.append([f'Train error ROM - HROM, Cluster {i}: {rom_hrom_train_error:.2E}, Modes: {modes}, Points: {len(mu)}'])
    
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    if Test_rom:
        KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('fom_rom_test_errors.npy')
        KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('rom_hrom_test_errors.npy')
        KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('test_rom_data.xlsx')
        KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('test_hrom_data.xlsx')

        if os.path.exists('predicted_indexes_list.npy'):
            predicted_indexes_list = np.int0(np.load('predicted_indexes_list.npy'))
        else:
            predicted_indexes_list = PEBL_Clustering_RBF_Prediction(    bisection_tolerance      = bisection_tolerance, 
                                                                        POD_tolerance            = POD_tolerance, 
                                                                        plot_clusters_prediction = True,
                                                                        plot_custering_info      = True,
                                                                        mu_train                 = mu_train,
                                                                        mu_train_not_scaled      = mu_train_not_scaled,
                                                                        mu_test                  = mu_test,
                                                                        mu_test_not_scaled       = mu_test_not_scaled)
            np.save('predicted_indexes_list', predicted_indexes_list)
    
        for i in range(np.int0(np.max(np.array(predicted_indexes_list))) + 1):

                # if i == cluster_number:   
                    
                    if not os.path.exists(f'Mu_by_clusters/mu_test_{i}.npy'):
                        new_mu_test = np.array([mu for mu, cluster in zip(mu_test, predicted_indexes_list) if cluster == i])
                        np.save(f'Mu_by_clusters/mu_test_{i}', new_mu_test)
                    else:
                        new_mu_test = np.load(f'Mu_by_clusters/mu_test_{i}.npy')

                    # number_of_mu_test = 10
                    # if number_of_mu_test < len(new_mu_test):
                    #     new_mu_test = np.array(random.sample(list(np.load(f'Mu_by_clusters/mu_test_{i}.npy')), number_of_mu_test))

                    # RBF PREDICTION
                    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                    mu = np.load(f'Mu_by_clusters/mu_train_{i}.npy')
                    mu_not_scaled = [mu_value for mu_value, mu_train_value in zip(mu_train_not_scaled, mu_train) if mu_train_value in mu]
                    new_mu_test_not_scaled = [mu_value for mu_value, mu_test_value in zip(mu_test_not_scaled, mu_test) if mu_test_value in new_mu_test]

                    RBF_prediction(mu_train            = mu,
                                   mu_train_not_scaled = mu_not_scaled, 
                                   mu_test             = new_mu_test, 
                                   mu_test_not_scaled  = new_mu_test_not_scaled)

                    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                    general_rom_manager_parameters = GetRomManagerParameters()
                    general_rom_manager_parameters['ROM']['rom_basis_output_folder'].SetString(f'rom_bases/rom_data_cluster_{i}')

                    project_parameters_name = "ProjectParameters.json"

                    rom_manager = RomManager(   project_parameters_name, 
                                                general_rom_manager_parameters, 
                                                CustomizeSimulation, 
                                                UpdateProjectParameters, 
                                                UpdateMaterialParametersFile )            
                    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

                    fom_rom_test_error, rom_hrom_test_error = rom_manager.Test(new_mu_test, store_all_snapshots = True)

                    modes = np.load(f'rom_bases/rom_data_cluster_{i}/RightBasisMatrix.npy').shape[1]
                    fom_rom_test_errors.append([f'Test  error FOM - ROM, cluster {i}: {fom_rom_test_error:.2E}, Modes: {modes}, Points: {len(new_mu_test)}'])
                    if rom_hrom_test_error != 0:
                        rom_hrom_test_errors.append([f'Test  error ROM - HROM, cluster {i}: {rom_hrom_test_error:.2E}, Modes: {modes}, Points: {len(new_mu_test)}'])


    print('::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
    PEBL_error_estimation(mu_train, mu_test)
    print('::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
    for i in range(len(fom_rom_train_errors)):
        print(fom_rom_train_errors[i])  
    print('::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
    for i in range(len(rom_hrom_train_errors)):
        print(rom_hrom_train_errors[i])  
    print('::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
    for i in range(len(fom_rom_test_errors)):
        print(fom_rom_test_errors[i])
    print('::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
    for i in range(len(rom_hrom_test_errors)):
        print(rom_hrom_test_errors[i])
    print('::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
    
    np.save('fom_rom_train_errors', fom_rom_train_errors)
    np.save('rom_hrom_train_errors', rom_hrom_train_errors)
    np.save('fom_rom_test_errors', fom_rom_test_errors)
    np.save('rom_hrom_test_errors', rom_hrom_test_errors)
