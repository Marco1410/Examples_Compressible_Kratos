import os
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
import numpy as np
import time as time
import KratosMultiphysics
import KratosMultiphysics.kratos_utilities
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp 
from local_rom_manager import LocalRomManager

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def CustomizeSimulation(cls, global_model, parameters, mesh_file_name, cluster, is_test=False):

    class CustomSimulation(cls):

        def __init__(self, model,project_parameters, custom_param = None):
            super().__init__(model,project_parameters)
            self.custom_param  = custom_param
        
        def Initialize(self):
            super().Initialize()
            self._GetSolver()._GetSolutionStrategy().SetKeepSystemConstantDuringIterations(False)

        def InitializeSolutionStep(self):
            super().InitializeSolutionStep()

            # data_name = f"Clustering/{mesh_file_name}/Centroids/{cluster}.npy"
            # u = np.load(data_name)

            # phi  = np.load(f"{mesh_file_name}/RomBases/rom_data_cluster_{cluster}/RightBasisMatrix.npy")
            # data = phi@(phi.T@u)

            # model_part = self.model["MainModelPart"]
            # for node in model_part.Nodes:
            #     offset = np.where(np.arange(1,model_part.NumberOfNodes()+1, dtype=int) == node.Id)[0][0]*2
            #     node.SetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL, data[offset+1])
            #     node.SetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL, data[offset])

        def FinalizeSolutionStep(self):
            super().FinalizeSolutionStep()
            if is_test:
                nametype = parameters["output_processes"]["vtk_output"][0]["Parameters"]["output_path"].GetString()
                simulation_name = nametype.removeprefix(f"{mesh_file_name}/Test/ROM_Results/")
                skin_data_filename = f"{mesh_file_name}/Test/ROM_Skin_Data/{simulation_name}.dat"
                capture_filename = f"{mesh_file_name}/Test/ROM_Captures/{simulation_name}.png"
            else:
                nametype = parameters["output_processes"]["vtk_output"][0]["Parameters"]["output_path"].GetString()
                simulation_name = nametype.removeprefix(f"{mesh_file_name}/ROM_Results/")
                skin_data_filename = f"{mesh_file_name}/ROM_Skin_Data/{simulation_name}.dat"
                capture_filename = f"{mesh_file_name}/ROM_Captures/{simulation_name}.png"
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
            x  = np.loadtxt(skin_data_filename, usecols=(0,))
            cp = np.loadtxt(skin_data_filename, usecols=(3,))
            fig = plt.plot(x, cp, 'xr', markersize = 2.0, label = simulation_name)
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

def UpdateProjectParameters(parameters, mu=None, mesh_file_name=None, is_test=False):
    angle_of_attack        = mu[0]
    mach_infinity          = mu[1]
    parameters["modelers"][0]["parameters"]["input_filename"].SetString(f'salome_mesh_files/{mesh_file_name}.med')
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(angle_of_attack)
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach_infinity)
    if is_test:
        parameters["output_processes"]["vtk_output"][0]["Parameters"]["output_path"].SetString(f'{mesh_file_name}/Test/ROM_Results/{angle_of_attack}, {mach_infinity}')
    else:
        parameters["output_processes"]["vtk_output"][0]["Parameters"]["output_path"].SetString(f'{mesh_file_name}/ROM_Results/{angle_of_attack}, {mach_infinity}')
    return parameters

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# save/load parameters
#
def save_mu_parameters(mu_train, name1, mu_test, name2):
    if len(mu_train) > 0:
        np.save(f'{name1}',mu_train)
    if len(mu_test) > 0:
        np.save(f'{name2}',mu_test)

def load_mu_parameters():
    if os.path.exists("mu_parameters/mu_train.npy") and os.path.exists("mu_parameters/mu_test.npy"):
        archivo = 'mu_parameters/mu_train.npy'
        mu_train = np.load(archivo)
        archivo = 'mu_parameters/mu_test.npy'
        mu_test = np.load(archivo)
    elif os.path.exists("mu_parameters/mu_train.npy"):
        archivo = 'mu_parameters/mu_train.npy'
        mu_train = np.load(archivo)
        mu_test = []
    elif os.path.exists("mu_parameters/mu_test.npy"):
        archivo = 'mu_parameters/mu_test.npy'
        mu_test = np.load(archivo)
        mu_train = []
    return np.array(mu_train), np.array(mu_test)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def GetRomManagerParameters():
    general_rom_manager_parameters = KratosMultiphysics.Parameters("""{
            "rom_stages_to_train" : [],      // ["ROM","HROM"]
            "rom_stages_to_test"  : ["ROM"],      // ["ROM","HROM"]
            "projection_strategy": "lspg",           
            "ROM":{
                "svd_truncation_tolerance": 1e-12,
                "model_part_name": "MainModelPart",
                "nodal_unknowns": ["VELOCITY_POTENTIAL","AUXILIARY_VELOCITY_POTENTIAL"], // Main unknowns. Snapshots are taken from these
                "rom_basis_output_format": "numpy",                                       // "json" "numpy"
                "rom_basis_output_name": "RomParameters",
                "rom_basis_output_folder": "rom_data",
                "snapshots_control_type": "step",                                        // "step", "time"
                "snapshots_interval": 1,
                "lspg_rom_bns_settings": {
                    "train_petrov_galerkin": false,
                    "basis_strategy": "reactions",                        // 'residuals', 'jacobian', 'reactions'
                    "include_phi": false,
                    "svd_truncation_tolerance": 1e-12,
                    "solving_technique": "normal_equations",              // 'normal_equations', 'qr_decomposition'
                    "monotonicity_preserving": false
                }
            }
        }""")
    return general_rom_manager_parameters
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

if __name__ == "__main__":

    update_bases = False

    with os.scandir('Clustering') as folders:
        folders = [dir.name for dir in folders if dir.is_dir() and '(copia)' not in dir.name]
        folders.sort()

    mu_train, mu_test = load_mu_parameters()
    mu_train = list(mu_train)
    mu_test  = list(mu_test)

    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    general_rom_manager_parameters = GetRomManagerParameters()
    project_parameters_name = "ProjectParametersFiles/RomProjectParameters.json"
    rom_manager = LocalRomManager(project_parameters_name,general_rom_manager_parameters,CustomizeSimulation,UpdateProjectParameters)
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    mesh_errors = []
    for mesh_file in folders:

        if not os.path.exists(f"{mesh_file}"):
            os.mkdir(f"{mesh_file}")

        mu_train_with_indexes = np.load(f'Clustering/{mesh_file}/mu_train_with_indexes.npy')
        predicted_indexes     = np.load(f'Clustering/{mesh_file}/predicted_indexes.npy')

        start_time = time.time()

        mu_train_errors = rom_manager.Fit(mu_train, mu_train_with_indexes, mesh_file, update_bases)

        exe_time = time.time() - start_time
        print('Executing ROM Manager took ' + str(round(exe_time, 2)) + ' sec')

        print(f'Mesh name | Approximation error FOM vs ROM')
        mesh_errors.append([mesh_file, rom_manager.PrintErrors()])
        print(mesh_errors)

        input("Pause")

        start_time = time.time()
        
        mu_test_errors  = rom_manager.Test(mu_test, mesh_file, predicted_indexes)

        exe_time = time.time() - start_time
        print('Executing ROM Manager took ' + str(round(exe_time, 2)) + ' sec')

        print(f'Mesh name | Approximation error FOM vs ROM')
        mesh_errors.append([mesh_file, rom_manager.PrintErrors()])
        print(mesh_errors)

        save_mu_parameters(mu_train_errors, f'{mesh_file}/mu_train_errors', mu_test_errors, f'{mesh_file}/mu_test_errors')