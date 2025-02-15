import os
import importlib
import numpy as np
import matplotlib.pyplot as plt
import KratosMultiphysics
import KratosMultiphysics.kratos_utilities
from KratosMultiphysics.RomApplication.rom_database import RomDatabase
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
from KratosMultiphysics.gid_output_process import GiDOutputProcess
from KratosMultiphysics.vtk_output_process import VtkOutputProcess


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Launch Fake Simulation
#
def LaunchFakeSimulation(data_set, mode, prefix):
    with open('ProjectParameters.json','r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    model = KratosMultiphysics.Model()
    parameters_copy = FakeProjectParameters(parameters.Clone(), mode, prefix)
    analysis_stage_class = _GetAnalysisStageClass(parameters_copy)
    simulation = FakeSimulation(analysis_stage_class, model, parameters_copy, data_set)
    simulation.Run()

    for process in simulation._GetListOfOutputProcesses():
            if isinstance(process, GiDOutputProcess):
                gid_output = process
    parameters_output = parameters_copy["output_processes"]["gid_output"][0]["Parameters"]['postprocess_parameters']
    gid_output = GiDOutputProcess(simulation.model['MainModelPart'],
                                    parameters_copy["output_processes"]["gid_output"][0]["Parameters"]['output_name'].GetString(),
                                    parameters_output)
    gid_output.ExecuteInitialize()
    gid_output.ExecuteBeforeSolutionLoop()
    gid_output.ExecuteInitializeSolutionStep()
    gid_output.PrintOutput()
    gid_output.ExecuteFinalizeSolutionStep()
    gid_output.ExecuteFinalize()

    for process in simulation._GetListOfOutputProcesses():
            if isinstance(process, VtkOutputProcess):
                vtk_output = process
    parameters_output = parameters_copy["output_processes"]["vtk_output"][0]["Parameters"]
    vtk_output = VtkOutputProcess(simulation.model, parameters_output)
    vtk_output.ExecuteInitialize()
    vtk_output.ExecuteBeforeSolutionLoop()
    vtk_output.ExecuteInitializeSolutionStep()
    vtk_output.PrintOutput()
    vtk_output.ExecuteFinalizeSolutionStep()
    vtk_output.ExecuteFinalize()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# save parameters
#
def FakeSimulation(cls, global_model, parameters, data_set):

    class CustomSimulation(cls):
        def __init__(self, model,project_parameters):
            super().__init__(model,project_parameters)

        def Run(self):
            self.Initialize()
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()
            self.Finalize()

        def Initialize(self):
            super().Initialize()
            model_part = self.model["MainModelPart"]
            for node in model_part.Nodes:
                offset = np.where(np.arange(1,model_part.NumberOfNodes()+1, dtype=int) == node.Id)[0][0]*2

                node.SetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL, data_set[offset])
                node.SetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL, data_set[offset+1])

        def InitializeSolutionStep(self):
            super().InitializeSolutionStep()

        def FinalizeSolutionStep(self):
            super().FinalizeSolutionStep()

    return CustomSimulation(global_model, parameters)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# save parameters
#
def FakeProjectParameters(parameters, mode, prefix):
    parameters["solver_settings"]["echo_level"].SetInt(0)
    parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString(f'Modes/{prefix}_{mode}')
    parameters["output_processes"]["vtk_output"][0]["Parameters"]["output_path"].SetString(f'Modes/{prefix}_{mode}')
    return parameters

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# save parameters
#
def _GetAnalysisStageClass(parameters):
    analysis_stage_module_name = parameters["analysis_stage"].GetString()
    analysis_stage_class_name = analysis_stage_module_name.split('.')[-1]
    analysis_stage_class_name = ''.join(x.title() for x in analysis_stage_class_name.split('_'))
    analysis_stage_module = importlib.import_module(analysis_stage_module_name)
    analysis_stage_class = getattr(analysis_stage_module, analysis_stage_class_name)
    return analysis_stage_class

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# load parameters
#
def load_mu_parameters(name):
    filename = f'{name}.npy'
    if os.path.exists(filename):
        mu_npy = np.load(filename)
        mu =  [mu.tolist() for mu in mu_npy]
    else:
        mu = []
    return mu

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def GetRomManagerParameters():
    general_rom_manager_parameters = KratosMultiphysics.Parameters("""{
            "rom_stages_to_train" : ["ROM"],            // ["ROM","HROM"]
            "rom_stages_to_test"  : ["ROM"],            // ["ROM","HROM"]
            "paralellism" : null,                       // null, TODO: add "compss"
            "projection_strategy": "galerkin",          // "lspg", "galerkin", "petrov_galerkin"
            "type_of_decoder" : "linear",               // "linear" "ann_enhanced",  TODO: add "quadratic"
            "assembling_strategy": "global",            // "global", "elemental"
            "save_gid_output": true,                    // false, true #if true, it must exits previously in the ProjectParameters.json
            "save_vtk_output": true,                   // false, true #if true, it must exits previously in the ProjectParameters.json
            "output_name": "id",                        // "id" , "mu"
            "store_nonconverged_fom_solutions": true,
            "ROM":{
                "svd_truncation_tolerance": 1e-12,
                "print_singular_values": true,
                "use_non_converged_sols" : false,
                "model_part_name": "MainModelPart",                         // This changes depending on the simulation: Structure, FluidModelPart, ThermalPart #TODO: Idenfity it automatically
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
                },
                "ann_enhanced_settings": {
                    "modes":[5,50],
                    "layers_size":[200,200],
                    "batch_size":2,
                    "epochs":800,
                    "NN_gradient_regularisation_weight": 0.0,
                    "lr_strategy":{
                        "scheduler": "sgdr",
                        "base_lr": 0.001,
                        "additional_params": [1e-4, 10, 400]
                    },
                    "training":{
                        "retrain_if_exists" : false  // If false only one model will be trained for each the mu_train and NN hyperparameters combination
                    },
                    "online":{
                        "model_number": 0   // out of the models existing for the same parameters, this is the model that will be lauched
                    }
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
                "include_nodal_neighbouring_elements_model_parts_list":[],
                "include_minimum_condition": false,
                "include_condition_parents": true
            }
        }""")

    return general_rom_manager_parameters

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


if __name__ == "__main__":

    modes_to_plot = [0, 1, 2, 3]
    strategies = ['galerkin','lspg']
    n = [0, 0]

    #Disable logs
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    
    for id, strategy in zip(n, strategies):

        prefix = f'{id}_{strategy}'
        mu_train = load_mu_parameters(f'Mu_history/{prefix}_mu_train')
        
        general_rom_manager_parameters = GetRomManagerParameters()
        general_rom_manager_parameters["projection_strategy"].SetString(strategy)
        data_base = RomDatabase(general_rom_manager_parameters, mu_names=None)

        is_sigma_in_database, hash_sigma = data_base.check_if_in_database("SingularValues_Solution", mu_train)
        if is_sigma_in_database:
            s = data_base.get_single_numpy_from_database(hash_sigma)

            valores_singulares = np.sort(s)[::-1]

            suma_acumulada = np.cumsum(valores_singulares)

            suma_acumulada_normalizada = suma_acumulada / suma_acumulada[-1]

            decaimiento = 1 - suma_acumulada_normalizada

            plt.figure(figsize=(10, 6))
            plt.semilogy(range(len(valores_singulares)), decaimiento, marker='x', markersize=5, linestyle='', color='b')
            plt.title('Singular values decay')
            plt.xlabel('Singular value index')
            plt.ylabel('Decay')
            plt.grid(True, which="both", ls="--")
            plt.savefig(f'Singular_values_decay_{prefix}.pdf', bbox_inches='tight')
            # plt.show()
            plt.close()

        is_basis_in_database, hash_basis = data_base.check_if_in_database("RightBasis", mu_train)
        if is_basis_in_database:
            
            phi = data_base.get_single_numpy_from_database(hash_basis)

            for mode in modes_to_plot:
                mode_values = phi[:,mode]
                LaunchFakeSimulation(mode_values, mode, prefix)
