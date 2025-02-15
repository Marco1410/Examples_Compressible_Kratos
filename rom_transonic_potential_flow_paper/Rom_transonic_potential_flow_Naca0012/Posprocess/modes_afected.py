import os
import numpy as np
import matplotlib.pyplot as plt
import KratosMultiphysics
import KratosMultiphysics.kratos_utilities
import KratosMultiphysics.CompressiblePotentialFlowApplication
from KratosMultiphysics.RomApplication.rom_database import RomDatabase

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

    # Desactivar logs innecesarios
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

    # Leer los IDs de los elementos desde el archivo .txt
    element_ids = np.loadtxt("selected_elements_list.txt", dtype=int)

    # Cargar el modelo y obtener la malla
    model = KratosMultiphysics.Model()
    model_part = model.CreateModelPart("MainModelPart")

    # Cargar la malla desde un archivo .mdpa o previamente generado
    KratosMultiphysics.ModelPartIO("Mesh/model_mesh_0").ReadModelPart(model_part)

    # Obtener los nodos de los elementos de la lista
    nodos_set = set()
    for elem_id in element_ids:
        element = model_part.GetElement(elem_id)
        for node in element.GetNodes():
            nodos_set.add(node.Id)

    # Convertir nodos únicos en un array ordenado
    nodos_lista = np.array(sorted(nodos_set))

    # Crear un diccionario de node_id -> índice en la lista de nodos ordenada
    node_id_to_index = {nid: idx for idx, nid in enumerate(nodos_lista)}

    # Obtener los índices de fila correspondientes a los nodos de los elementos seleccionados
    row_indexes = np.array([node_id_to_index[nid] for nid in nodos_lista])

    strategies = ['galerkin']
    modes_to_plot = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
    n = [4]

    constan = 0.01

    for id, strategy in zip(n, strategies):

        prefix = f'{strategy}'
        # mu_train = load_mu_parameters(f'Mu_history/{prefix}_mu_train')
        mu_train = load_mu_parameters('Mu_history/mu_train')
        
        general_rom_manager_parameters = GetRomManagerParameters()
        general_rom_manager_parameters["projection_strategy"].SetString(strategy)
        data_base = RomDatabase(general_rom_manager_parameters, mu_names=None)

        is_basis_in_database, hash_basis = data_base.check_if_in_database("RightBasis", mu_train)
        
        if is_basis_in_database:
            phi = data_base.get_single_numpy_from_database(hash_basis)

            # Generar un array de índices de todas las filas
            all_rows = np.arange(phi.shape[0])

            # Filtrar las filas que NO están en `row_indexes`
            rows_to_modify = np.setdiff1d(all_rows, row_indexes)  # Complemento de `row_indexes`

            # phi[row_indexes,:] *= constan
            phi[rows_to_modify,:] *= constan

            # Graficar los modos (columnas de phi)
            plt.figure(figsize=(12, 6))

            for mode in modes_to_plot:
                plt.plot(phi[:, mode], label=f'Modo {mode+1}')

            plt.title("Primeros Modos del Flujo Transónico (POD)")
            plt.legend()
            # plt.savefig('First_POD_modes.png')
            plt.show()
            plt.close()