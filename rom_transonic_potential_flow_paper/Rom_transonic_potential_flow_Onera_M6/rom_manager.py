import os
import time
import openpyxl
import random
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
from plot import Plot_Cps # type: ignore
from ClearAll import Clear # type: ignore
from scipy.stats import qmc  
import KratosMultiphysics
import KratosMultiphysics.kratos_utilities
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
from KratosMultiphysics.RomApplication.rom_manager import RomManager
from KratosMultiphysics.RomApplication.calculate_rom_basis_output_process import CalculateRomBasisOutputProcess


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# get multiple parameters by Halton or LatinHypercube methods
#
def get_multiple_parameters(number_train_values=0, number_test_values=0, angle=[], mach=[], mu_validation=[None], method='Halton'):
    if method == 'Halton':
        sampler = qmc.Halton(d=2)
    elif method == 'LatinHypercube':
        sampler = qmc.LatinHypercube(d=2)
    mu_train = []; mu_test = []
    mu_train_not_scaled = []; mu_test_not_scaled = []; mu_validation_not_scaled = []
    if number_train_values > 0:
        sample = sampler.random(number_train_values)
        values = qmc.scale(sample, [angle[0],mach[0]], [angle[1],mach[1]])
        for i in range(number_train_values):
            #Angle of attack , Mach infinit
            mu_train.append([values[i,0], values[i,1]])
            mu_train_not_scaled.append([sample[i,0], sample[i,1]])
        np.save(f'mu_train', mu_train)
        np.save(f'mu_train_not_scaled', mu_train_not_scaled)
    if number_test_values > 0:
        sample = sampler.random(number_test_values)
        values = qmc.scale(sample, [angle[0],mach[0]], [angle[1],mach[1]])
        for i in range(number_test_values):
            #Angle of attack , Mach infinit
            mu_test.append([values[i,0], values[i,1]])
            mu_test_not_scaled.append([sample[i,0], sample[i,1]])
        np.save(f'mu_test', mu_test)
        np.save(f'mu_test_not_scaled', mu_test_not_scaled)
    if len(mu_validation)>0:
        mu_validation_not_scaled = get_not_scale_parameters(mu_validation, angle, mach)
        np.save(f'mu_validation', mu_validation)
        np.save(f'mu_validation_not_scaled', mu_validation_not_scaled)
    
    plot_mu_values(mu_train, mu_test, mu_validation, 'MuValues')
    plot_mu_values(mu_train_not_scaled, mu_test_not_scaled, mu_validation_not_scaled, 'MuValuesNotScaled')
    
    return mu_train, mu_test, mu_train_not_scaled, mu_test_not_scaled, mu_validation_not_scaled

def get_not_scale_parameters(mu=[None], angle=[], mach=[]):
    mu_not_scaled = []
    if len(mu) > 0:
        sample = np.array(mu)
        values = qmc.scale(sample, [angle[0],mach[0]], [angle[1],mach[1]], reverse=True)
        for i in range(len(mu)):
            #Angle of attack , Mach infinit
            mu_not_scaled.append([values[i,0], values[i,1]])    
    return mu_not_scaled


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# plot parameters
#
def plot_mu_values(mu_train, mu_test, mu_validation, name):
    if len(mu_train) > 0: plt.plot(np.array(mu_train)[:,1], np.array(mu_train)[:,0], 'bs', label="Train Values")
    if len(mu_test ) > 0: plt.plot(np.array( mu_test)[:,1], np.array( mu_test)[:,0], 'ro', label="Test Values")
    if len(mu_validation ) > 0: plt.plot(np.array( mu_validation)[:,1], np.array( mu_validation)[:,0], 'g*', label="Validation Values")
    plt.title('Mu Values')
    plt.ylabel('Alpha')
    plt.xlabel('Mach')
    plt.grid(True)
    plt.legend(bbox_to_anchor=(0.78, 1.05, 1.0, 0.1), loc='upper left', borderaxespad=0.)
    plt.savefig(f"{name}.png")
    plt.close('all')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# load parameters
#
def load_mu_parameters():
    if os.path.exists("mu_train.npy"):
        mu_train = np.load('mu_train.npy')
        mu_train_not_scaled = np.load('mu_train_not_scaled.npy')
        mu_train =  [mu.tolist() for mu in mu_train]
        mu_train_not_scaled =  [mu.tolist() for mu in mu_train_not_scaled]
    else:
        mu_train = []
        mu_train_not_scaled = []
    if os.path.exists("mu_test.npy"):
        mu_test = np.load('mu_test.npy')
        mu_test_not_scaled = np.load('mu_test_not_scaled.npy')
        mu_test =  [mu.tolist() for mu in mu_test]
        mu_test_not_scaled =  [mu.tolist() for mu in mu_test_not_scaled]
    else:
        mu_test = []
        mu_test_not_scaled = []
    if os.path.exists("mu_validation.npy"):
        mu_validation = np.load('mu_validation.npy')
        mu_validation_not_scaled = np.load('mu_validation_not_scaled.npy')
        mu_validation =  [mu.tolist() for mu in mu_validation]
        mu_validation_not_scaled =  [mu.tolist() for mu in mu_validation_not_scaled]
    else:
        mu_validation = []
        mu_validation_not_scaled = []
    return mu_train, mu_test, mu_train_not_scaled, mu_test_not_scaled, mu_validation, mu_validation_not_scaled

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# RBF prediction and error
#
def RBF_prediction( mu_train                 = [None],
                    mu_test                  = [None],
                    mu_train_not_scaled      = [None],
                    mu_test_not_scaled       = [None]):
    from ezyrb import ReducedOrderModel as ROM
    from ezyrb import RBF, POD, Database

    #### CLUSTERING DATA
    #################################################################################################
    parameters = np.array(mu_train_not_scaled)

    snapshots = []
    for mu in mu_train:
        file = f'{mu[0]}, {mu[1]}.dat'
        snapshots.append(np.array(np.loadtxt(f'FOM_Skin_Data/{file}', usecols=(3,))).reshape(-1,1))
    snapshots = np.block(snapshots)

    #### RBF TRAINING
    #################################################################################################
    db = Database(parameters, snapshots.T)
    pod = POD()
    rbf = RBF()
    rom = ROM(db, pod, rbf).fit()

    if len(mu_test) > 0:
        #### PREDICTION OF TEST
        #################################################################################################
        interpolated_solutions_list = [rom.predict([element]).snapshots_matrix for element in mu_test_not_scaled]

        for i, solution in enumerate(interpolated_solutions_list):
            np.save(f"RBF_Skin_Data/{mu_test[i][0]}, {mu_test[i][1]}", solution.T)

def RBF_error_estimation(mu_train, mu_test):
    if len(mu_train) > 0:
        approximation_error = 0.0
        FOM_model = []; RBF_model = []
        for mu in mu_train:
            FOM_model.append(np.array(np.loadtxt(f'FOM_Skin_Data/{mu[0]}, {mu[1]}.dat', usecols=(3,))).reshape(-1,1))
            RBF_model.append(np.array(np.load(f"RBF_Skin_Data/{mu[0]}, {mu[1]}.npy")).reshape(-1,1))
        FOM_model = np.block(FOM_model)
        RBF_model = np.block(RBF_model)
        training_approximation_error = np.linalg.norm(FOM_model - RBF_model)/np.linalg.norm(FOM_model)
        print(f'RBF training approximation error: {training_approximation_error:.2E}')

    if len(mu_test)>0 and os.path.exists(f'RBF_Skin_Data/{mu_test[0][0]}, {mu_test[0][1]}.npy'):
        approximation_error = 0.0
        FOM_model = []; RBF_model_interpolation = []
        for mu in mu_test:
            FOM_model.append(np.array(np.loadtxt(f'FOM_Skin_Data/{mu[0]}, {mu[1]}.dat', usecols=(3,))).reshape(-1,1))
            RBF_model_interpolation.append(np.array(np.load(f"RBF_Skin_Data/{mu[0]}, {mu[1]}.npy")).reshape(-1,1))
        FOM_model = np.block(FOM_model)
        RBF_model_interpolation = np.block(RBF_model_interpolation)
        approximation_error = np.linalg.norm(FOM_model - RBF_model_interpolation)/np.linalg.norm(FOM_model)
        print(f'RBF interpolation approximation error: {approximation_error:.2E}')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


def CustomizeSimulation(cls, global_model, parameters, mu):

    class CustomSimulation(cls):

        def __init__(self, model,project_parameters):
            super().__init__(model,project_parameters)
        
        def Initialize(self):
            # if self._GetSimulationName() == "::[ROM Simulation]:: ":
            #     parameters["solver_settings"]["solving_strategy_settings"]["advanced_settings"]["min_alpha"].SetDouble(0.05)
            #     parameters["solver_settings"]["solving_strategy_settings"]["advanced_settings"]["line_search_tolerance"].SetDouble(0.25)
            super().Initialize()

            # self._GetSolver()._GetSolutionStrategy().SetKeepSystemConstantDuringIterations(True)
                
            print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
            print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
            print(f"{self._GetSimulationName()}")
            angle = parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].GetDouble()
            mach  = parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].GetDouble()
            case_name = f'{angle}, {mach}'
            print(f"{case_name}")
            if parameters["output_processes"].Has("gid_output"):
                print(f"{parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].GetString().removeprefix('Results/')}")
            print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
            print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")

        def Run(self):

            angle = parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].GetDouble()
            mach  = parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].GetDouble()

            case_name = f'{angle}, {mach}'

            info_steps_list = []
            error = 0
            modes = 0

            if self._GetSimulationName() == "Analysis": # FOM

                start_time = time.time()
                self.Initialize()
                self.RunSolutionLoop()
                self.Finalize()

                exe_time = time.time() - start_time

                if parameters["output_processes"].Has("gid_output"):

                    simulation_name = parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].GetString().removeprefix('Results/')

                    for process in self._GetListOfOutputProcesses():
                            if isinstance(process, CalculateRomBasisOutputProcess):
                                BasisOutputProcess = process

                    if 'Fit' in simulation_name:
                        case_type = 'train_fom'
                    elif 'Test' in simulation_name:
                        case_type = 'test_fom' 
                        modes = np.load('rom_data/RightBasisMatrix.npy').shape[1]
                    elif 'Run' in simulation_name:
                        case_type = 'run_fom' 
                    skin_data_filename = f"FOM_Skin_Data/{case_name}.dat"
                    fom = BasisOutputProcess._GetSnapshotsMatrix()
                    np.save(f'FOM_Snapshots/{case_name}',fom)

                    fout = open(skin_data_filename,'w')
                    modelpart = self.model["MainModelPart.Wing"]
                    for node in modelpart.Nodes:
                        x = node.X ; y = node.Y ; z = node.Z
                        cp = node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
                        fout.write("%s %s %s %s\n" %(x,y,z,cp))
                    fout.close()

                    modelpart = self.model["MainModelPart"]
                    fout = open("trailing_edge_element_id.txt",'a')
                    for elem in modelpart.Elements:
                        if elem.GetValue(CPFApp.TRAILING_EDGE):
                            fout.write("%s\n" %(elem.Id))
                        if elem.GetValue(CPFApp.KUTTA):
                            fout.write("%s\n" %(elem.Id))
                    fout.close()

                    fout = open("upwind_elements_list.txt",'a')
                    for elem in modelpart.Elements:
                        if elem.GetValue(CPFApp.ID_UPWIND_ELEMENT) != 0.0:
                            fout.write("%s\n" %(elem.Id))
                            fout.write("%s\n" %(elem.GetValue(CPFApp.ID_UPWIND_ELEMENT)))
                    fout.close()
                
                    info_steps_list.append([case_type,
                                            angle,
                                            mach, 
                                            self.model["MainModelPart"].ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER],
                                            self.model["MainModelPart"].ProcessInfo[KratosMultiphysics.RESIDUAL_NORM],
                                            error,
                                            modes,
                                            round(exe_time, 2),
                                            self.model["MainModelPart"].NumberOfNodes(),
                                            self.model["MainModelPart"].ProcessInfo[CPFApp.LIFT_COEFFICIENT],
                                            self.model["MainModelPart"].ProcessInfo[KratosMultiphysics.DRAG_COEFFICIENT]])
                    
                    if os.path.exists(f'FOM_data.xlsx'):
                        wb = openpyxl.load_workbook(f'FOM_data.xlsx')
                        hoja = wb.active
                        for item in info_steps_list:
                            hoja.append(item)
                        wb.save(f'FOM_data.xlsx')
                    else:
                        wb = openpyxl.Workbook()
                        hoja = wb.active
                        hoja.append(('Case name', 'Angle [º]', 'Mach', 'NL iterations', 'Residual norm', 'Approximation error [%]', 'Modes', 'Time [sec]', 'Nodes', 'Cl', 'Cd'))
                        for item in info_steps_list:
                            hoja.append(item)
                        wb.save(f'FOM_data.xlsx')

            elif self._GetSimulationName() == "::[ROM Simulation]:: ": # Global ROM

                start_time = time.time()
                self.Initialize()
                self.RunSolutionLoop()
                self.Finalize()                
                exe_time = time.time() - start_time

                if parameters["output_processes"].Has("gid_output"):

                    simulation_name = parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].GetString().removeprefix('Results/')

                    for process in self._GetListOfOutputProcesses():
                            if isinstance(process, CalculateRomBasisOutputProcess):
                                BasisOutputProcess = process

                    if 'HROM' in simulation_name:                      # HROM
                        if 'Fit' in simulation_name:
                            case_type = 'train_'
                        elif 'Test' in simulation_name:
                            case_type = 'test_'
                        elif 'Run' in simulation_name:
                            case_type = 'run_' 
                            
                        hrom = BasisOutputProcess._GetSnapshotsMatrix()
                        fom = np.load(f'FOM_Snapshots/{case_name}.npy')

                        if (len(fom) != len(hrom) or 'HHROM' in simulation_name):
                            q_matrix = []
                            main_model_part = self.model["MainModelPart"]
                            q_matrix.append(np.array(main_model_part.GetValue(KratosMultiphysics.RomApplication.ROM_SOLUTION_INCREMENT)).reshape(-1,1))
                            q_matrix = np.block(q_matrix)
                            phi = np.load(f'rom_data/RightBasisMatrix.npy')
                            if (q_matrix.shape[0] == 1): q_matrix = q_matrix.T
                            hrom = phi @ q_matrix
                            np.save(f'HHROM_Snapshots/{case_name}', hrom)
                            case_type = case_type + 'hhrom'
                        else:
                            np.save(f'HROM_Snapshots/{case_name}', hrom)
                            case_type = case_type + 'hrom'
                            skin_data_filename = f"HROM_Skin_Data/{case_name}.dat"
                            fout = open(skin_data_filename,'w')
                            modelpart = self.model["MainModelPart.Wing"]
                            for node in modelpart.Nodes:
                                x = node.X ; y = node.Y ; z = node.Z
                                cp = node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
                                fout.write("%s %s %s %s\n" %(x,y,z,cp))
                            fout.close()

                        error = np.linalg.norm(fom-hrom)/np.linalg.norm(fom)
                        modes = np.load('rom_data/RightBasisMatrix.npy').shape[1]
                    
                        info_steps_list.append([case_type,
                                                angle,
                                                mach, 
                                                self.model["MainModelPart"].ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER],
                                                self.model["MainModelPart"].ProcessInfo[KratosMultiphysics.RESIDUAL_NORM],
                                                error,
                                                modes,
                                                round(exe_time, 2),
                                                self.model["MainModelPart"].NumberOfNodes(),
                                                self.model["MainModelPart"].ProcessInfo[CPFApp.LIFT_COEFFICIENT],
                                                self.model["MainModelPart"].ProcessInfo[KratosMultiphysics.DRAG_COEFFICIENT]])
                        
                        if os.path.exists(f'case_data.xlsx'):
                            wb = openpyxl.load_workbook(f'case_data.xlsx')
                            hoja = wb.active
                            for item in info_steps_list:
                                hoja.append(item)
                            wb.save(f'case_data.xlsx')
                        else:
                            wb = openpyxl.Workbook()
                            hoja = wb.active
                            hoja.append(('Case name', 'Angle [º]', 'Mach', 'NL iterations', 'Residual norm', 'Approximation error [%]', 'Modes', 'Time [sec]', 'Nodes', 'Cl', 'Cd'))
                            for item in info_steps_list:
                                hoja.append(item)
                            wb.save(f'case_data.xlsx')

                    elif 'ROM' in simulation_name:                    # ROM
                        if 'Fit' in simulation_name:
                            case_type = 'train_rom'
                        elif 'Test' in simulation_name:
                            case_type = 'test_rom' 
                        elif 'Run' in simulation_name:
                            case_type = 'run_rom' 
                        modes = np.load('rom_data/RightBasisMatrix.npy').shape[1]
                        skin_data_filename = f"ROM_Skin_Data/{case_name}.dat"

                        rom = BasisOutputProcess._GetSnapshotsMatrix()
                        np.save(f'ROM_Snapshots/{case_name}',rom)
                        fom = np.load(f'FOM_Snapshots/{case_name}.npy')
                        error = np.linalg.norm(fom-rom)/np.linalg.norm(fom)

                        fout = open(skin_data_filename,'w')
                        modelpart = self.model["MainModelPart.Wing"]
                        for node in modelpart.Nodes:
                            x = node.X ; y = node.Y ; z = node.Z
                            cp = node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
                            fout.write("%s %s %s %s\n" %(x,y,z,cp))
                        fout.close()
                    
                        info_steps_list.append([case_type,
                                                angle,
                                                mach, 
                                                self.model["MainModelPart"].ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER],
                                                self.model["MainModelPart"].ProcessInfo[KratosMultiphysics.RESIDUAL_NORM],
                                                error,
                                                modes,
                                                round(exe_time, 2),
                                                self.model["MainModelPart"].NumberOfNodes(),
                                                self.model["MainModelPart"].ProcessInfo[CPFApp.LIFT_COEFFICIENT],
                                                self.model["MainModelPart"].ProcessInfo[KratosMultiphysics.DRAG_COEFFICIENT]])
                        
                        if os.path.exists(f'case_data.xlsx'):
                            wb = openpyxl.load_workbook(f'case_data.xlsx')
                            hoja = wb.active
                            for item in info_steps_list:
                                hoja.append(item)
                            wb.save(f'case_data.xlsx')
                        else:
                            wb = openpyxl.Workbook()
                            hoja = wb.active
                            hoja.append(('Case name','Angle [º]', 'Mach', 'NL iterations', 'Residual norm', 'Approximation error [%]', 'Modes', 'Time [sec]', 'Nodes', 'Cl', 'Cd'))
                            for item in info_steps_list:
                                hoja.append(item)
                            wb.save(f'case_data.xlsx')

    return CustomSimulation(global_model, parameters)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def UpdateProjectParameters(parameters, mu=None):
    angle_of_attack = mu[0]
    mach_infinity   = mu[1]
    wake_normal     = [-np.sin(angle_of_attack*np.pi/180),0.0,np.cos(angle_of_attack*np.pi/180)]
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(np.double(angle_of_attack))
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(np.double(mach_infinity))
    parameters["processes"]["boundary_conditions_process_list"][1]["Parameters"]["wake_process_cpp_parameters"]["wake_normal"].SetVector(wake_normal)

    return parameters

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def UpdateMaterialParametersFile(material_parametrs_file_name, mu):
    pass

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



def GetRomManagerParameters():
    general_rom_manager_parameters = KratosMultiphysics.Parameters("""{
            "rom_stages_to_train" : ["ROM","HROM"],            // ["ROM","HROM"]
            "rom_stages_to_test"  : ["ROM","HROM"],            // ["ROM","HROM"]
            "paralellism" : null,                       // null, TODO: add "compss"
            "projection_strategy": "galerkin",          // "lspg", "galerkin", "petrov_galerkin"
            "type_of_decoder" : "linear",               // "linear" "ann_enhanced",  TODO: add "quadratic"
            "assembling_strategy": "global",            // "global", "elemental"
            "save_gid_output": true,                    // false, true #if true, it must exits previously in the ProjectParameters.json
            "save_vtk_output": true,                   // false, true #if true, it must exits previously in the ProjectParameters.json
            "output_name": "id",                        // "id" , "mu"
            "store_nonconverged_fom_solutions": true,
            "ROM":{
                "analysis_stage" : "KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis",
                "svd_truncation_tolerance": 1e-6,
                "print_singular_values": true,
                "use_non_converged_sols" : false,
                "model_part_name": "MainModelPart",                            // This changes depending on the simulation: Structure, FluidModelPart, ThermalPart #TODO: Idenfity it automatically
                "nodal_unknowns": ["VELOCITY_POTENTIAL", "AUXILIARY_VELOCITY_POTENTIAL"],     // Main unknowns. Snapshots are taken from these
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
                    "svd_truncation_tolerance": 0,
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
                "create_hrom_visualization_model_part" : true,
                "echo_level" : 1,                                       
                "hrom_format": "numpy",
                "initial_candidate_elements_model_part_list" : [],
                "initial_candidate_conditions_model_part_list" : [],
                "include_elements_model_parts_list": [],
                "include_conditions_model_parts_list": [],
                "include_nodal_neighbouring_elements_model_parts_list":[],
                "include_minimum_condition": false,
                "include_condition_parents": true
            }
        }""")

    return general_rom_manager_parameters


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


if __name__ == "__main__":

    folder_names = ["FOM_Snapshots"  , "FOM_Skin_Data", 
                    "ROM_Snapshots"  , "ROM_Skin_Data", 
                    "HROM_Snapshots" , "HROM_Skin_Data",
                    "HHROM_Snapshots", "HHROM_Skin_Data", 
                    "RBF_Snapshots"  , "RBF_Skin_Data", 
                    "Train_Captures" , "Test_Captures", "Validation"]
    for name in folder_names:
        if not os.path.exists(name):
            os.mkdir(name)

    ###############################
    # PARAMETERS SETTINGS
    update_parameters  = True
    number_of_mu_train = 3
    number_of_mu_test  = 1
    re_dim_mu_train    = 0
    re_dim_mu_test     = 0
    delete_new_mu_val  = True
    angle_range        = [ 2.50, 3.25]
    mach_range         = [ 0.80, 0.85]
    ###############################

    mu_validation = []
    mu_validation_not_scaled = []
    mu_validation.append([3.06,0.839])
    
    if update_parameters:
        (mu_train, mu_test,
         mu_train_not_scaled, mu_test_not_scaled,
         mu_validation_not_scaled) = get_multiple_parameters(number_train_values = number_of_mu_train,
                                                            number_test_values   = number_of_mu_test , 
                                                            angle                = angle_range       , 
                                                            mach                 = mach_range        , 
                                                            mu_validation        = mu_validation     ,
                                                            method               = 'LatinHypercube'  )
        KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('upwind_elements_list.txt')
        KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('trailing_edge_element_id.txt')
        KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_train_new.npy')
        KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_train_not_scaled_new.npy')
        KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_test_new.npy')
        KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_test_not_scaled_new.npy')
    else:
        KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('case_data.xlsx')
        KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Train_Captures')
        KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Test_Captures')
        KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Validation')
        if delete_new_mu_val:
            KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_train_new.npy')
            KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_train_not_scaled_new.npy')
            KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_test_new.npy')
            KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_test_not_scaled_new.npy')
        for name in folder_names:
            if not os.path.exists(name):
                os.mkdir(name)
        if os.path.exists('mu_train_new.npy'):
            mu_train = np.load('mu_train_new.npy')
            mu_train =  [mu.tolist() for mu in mu_train]
            mu_train_not_scaled = np.load('mu_train_not_scaled_new.npy')
            mu_train_not_scaled =  [mu.tolist() for mu in mu_train_not_scaled]
            mu_test = np.load('mu_test_new.npy')
            mu_test =  [mu.tolist() for mu in mu_test]
            mu_test_not_scaled = np.load('mu_test_not_scaled_new.npy')
            mu_test_not_scaled =  [mu.tolist() for mu in mu_test_not_scaled]
        else:
            mu_train_list, mu_test_list, mu_train_not_scaled_list, mu_test_not_scaled_list, mu_validation, mu_validation_not_scaled = load_mu_parameters()
            mu_train = list(random.sample(mu_train_list, re_dim_mu_train))
            mu_train_not_scaled = [mu_value for mu_value, mu_train_value in zip(mu_train_not_scaled_list, mu_train) if mu_train_value in mu_train_list]
            mu_test = list(random.sample(mu_test_list, re_dim_mu_test))
            mu_test_not_scaled = [mu_value for mu_value, mu_test_value in zip(mu_test_not_scaled_list, mu_test) if mu_test_value in mu_test_list]
            np.save('mu_train_new',mu_train)
            np.save('mu_train_not_scaled_new',mu_train_not_scaled)
            np.save('mu_test_new',mu_test)
            np.save('mu_test_not_scaled_new',mu_test_not_scaled)
        plot_mu_values(mu_train, mu_test, mu_validation, 'MuValues')
        plot_mu_values(mu_train_not_scaled, mu_test_not_scaled, mu_validation_not_scaled, 'MuValuesNotScaled')

    general_rom_manager_parameters = GetRomManagerParameters()
    project_parameters_name = "ProjectParameters.json"

    rom_manager = RomManager(project_parameters_name,general_rom_manager_parameters,
                             CustomizeSimulation,UpdateProjectParameters,UpdateMaterialParametersFile,
                             relaunch_FOM=False, relaunch_ROM=False, relaunch_HROM=False)

    rom_manager.Fit(mu_train)

    rom_manager.Test(mu_test)

    rom_manager.RunFOM(mu_run=mu_validation)
    training_stages = rom_manager.general_rom_manager_parameters["rom_stages_to_train"].GetStringArray()
    if any(item == "ROM" for item in training_stages):
        rom_manager.RunROM(mu_run=mu_validation)
    if any(item == "HROM" for item in training_stages):
        rom_manager.RunHROM(mu_run=mu_validation,use_full_model_part=True)
    if any(item == "HHROM" for item in training_stages):
        rom_manager.RunHHROM(mu_run=mu_validation)

    if number_of_mu_train >= 3:
        RBF_prediction(mu_train = mu_train, mu_train_not_scaled = mu_train_not_scaled, 
                    mu_test = mu_train + mu_test, mu_test_not_scaled  = mu_train_not_scaled + mu_test_not_scaled)
        RBF_prediction(mu_train = mu_train, mu_train_not_scaled = mu_train_not_scaled, 
                    mu_test = mu_validation, mu_test_not_scaled  = mu_validation_not_scaled)

    Plot_Cps(mu_train, 'Train_Captures')
    Plot_Cps(mu_test, 'Test_Captures')
    Plot_Cps(mu_validation, 'Validation')

    rom_manager.PrintErrors(show_q_errors=False)

    if number_of_mu_train >= 3:
        print('::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
        RBF_error_estimation(mu_train, mu_test)
        print('::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
        print('Validation Error::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
        print('::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
        RBF_error_estimation(mu_train, mu_validation)
        print('::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::')