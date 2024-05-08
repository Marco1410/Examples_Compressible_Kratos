import os
import random
import openpyxl
import importlib
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
import numpy as np
import time as time
from ezyrb import Database
from ezyrb import ReducedOrderModel as ROM
from ezyrb import RBF, POD
import KratosMultiphysics
import KratosMultiphysics.kratos_utilities
from local_rom_manager import LocalRomManager
from KratosMultiphysics.RomApplication.calculate_rom_basis_output_process import CalculateRomBasisOutputProcess

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Custom ROM - HROM
#
def CustomizeSimulation(cls, global_model, parameters, mesh_name, cluster, is_test=False):

    class CustomSimulation(cls):

        def __init__(self, model,project_parameters, custom_param = None):
            super().__init__(model,project_parameters)
            self.custom_param  = custom_param
        
        def Initialize(self):
            super().Initialize()

        def InitializeSolutionStep(self):
            super().InitializeSolutionStep()

        def FinalizeSolutionStep(self):
            super().FinalizeSolutionStep()
            if is_test:
                nametype = parameters["output_processes"]["vtk_output"][0]["Parameters"]["output_path"].GetString()
                simulation_name = nametype.removeprefix(f"{mesh_name}/Test/ROM_Results/")
                skin_data_filename = f"{mesh_name}/Test/ROM_Skin_Data/{simulation_name}.dat"
                capture_filename = f"{mesh_name}/Test/ROM_Captures/{simulation_name}.png"
            else:
                nametype = parameters["output_processes"]["vtk_output"][0]["Parameters"]["output_path"].GetString()
                simulation_name = nametype.removeprefix(f"{mesh_name}/ROM_Results/")
                skin_data_filename = f"{mesh_name}/ROM_Skin_Data/{simulation_name}.dat"
                capture_filename = f"{mesh_name}/ROM_Captures/{simulation_name}.png"
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
        
        def Finalize(self):
            super().Finalize()

    return CustomSimulation(global_model, parameters)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def UpdateProjectParameters(parameters, mu=None, mesh_name=None, is_test=False):
    angle_of_attack        = mu[0]
    mach_infinity          = mu[1]
    parameters["modelers"][0]["parameters"]["input_filename"].SetString(f'{mesh_name}.med')
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(angle_of_attack)
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach_infinity)
    if is_test:
        parameters["output_processes"]["vtk_output"][0]["Parameters"]["output_path"].SetString(f'{mesh_name}/Test/ROM_Results/{angle_of_attack}, {mach_infinity}')
    else:
        parameters["output_processes"]["vtk_output"][0]["Parameters"]["output_path"].SetString(f'{mesh_name}/ROM_Results/{angle_of_attack}, {mach_infinity}')
    return parameters

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Custom FOM
#
def LaunchFOMSimulations(mu_test, mesh_file_name):
    if not os.path.exists(f"DataBase/Snapshots"):
        os.mkdir(f"DataBase/Snapshots")
    if not os.path.exists(f"DataBase/Snapshots_not_converged_steps"):
        os.mkdir(f"DataBase/Snapshots_not_converged_steps")
    if not os.path.exists(f"DataBase/FOM_Results"):
        os.mkdir(f"DataBase/FOM_Results")
    if not os.path.exists(f"DataBase/FOM_Skin_Data"):
        os.mkdir(f"DataBase/FOM_Skin_Data")
    if not os.path.exists(f"DataBase/FOM_Captures"):
        os.mkdir(f"DataBase/FOM_Captures")
    with open('FOMProjectParameters.json','r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    for Id, mu in enumerate(mu_test):
        fom_info_steps_list = []
        model = KratosMultiphysics.Model()
        parameters_copy = UpdateProjectParametersFOM(parameters.Clone(), mu, mesh_file_name)
        analysis_stage_class = _GetAnalysisStageClass(parameters_copy)
        simulation = CustomizeSimulationFOM(analysis_stage_class, model, parameters_copy, mesh_file_name)
        simulation.Run()
        model_part = simulation.model[parameters_copy["solver_settings"]["model_part_name"].GetString()].GetRootModelPart()
        fom_info_steps_list.append([Id,
                                    mu[0],
                                    mu[1], 
                                    model_part.ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER],
                                    model_part.ProcessInfo[KratosMultiphysics.RESIDUAL_NORM],
                                    ])
        for process in simulation._GetListOfOutputProcesses():
            if isinstance(process, CalculateRomBasisOutputProcess):
                BasisOutputProcess = process
        np.save(f'DataBase/Snapshots/{mu[0]}, {mu[1]}',BasisOutputProcess._GetSnapshotsMatrix()) 
        np.save(f'DataBase/Snapshots_not_converged_steps/{mu[0]}, {mu[1]}', np.array(simulation._GetSolver()._GetSolutionStrategy().GetIntermediateSolutionsMatrix()))

        if os.path.exists(f"DataBase/fom_data.xlsx"):
            wb = openpyxl.load_workbook(f"DataBase/fom_data.xlsx")
            hoja = wb.active
            for item in fom_info_steps_list:
                hoja.append(item)
            wb.save(f'DataBase/fom_data.xlsx')
        else:
            wb = openpyxl.Workbook()
            hoja = wb.active
            hoja.append(('Id', 'Angle', 'Mach', 'NL iterations', 'Residual norm'))
            for item in fom_info_steps_list:
                hoja.append(item)
            wb.save(f'DataBase/fom_data.xlsx')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def CustomizeSimulationFOM(cls, global_model, parameters, mesh_file_name):

    class CustomSimulation(cls):
        def __init__(self, model,project_parameters, custom_param = None):
            super().__init__(model,project_parameters)
            self.custom_param  = custom_param

        def Initialize(self):
            super().Initialize()

        def InitializeSolutionStep(self):
            super().InitializeSolutionStep()

        def FinalizeSolutionStep(self):
            super().FinalizeSolutionStep()
            nametype = parameters["output_processes"]["vtk_output"][0]["Parameters"]["output_path"].GetString()
            simulation_name = nametype.removeprefix(f"DataBase/FOM_Results/")
            skin_data_filename = f"DataBase/FOM_Skin_Data/{simulation_name}.dat"
            capture_filename = f"DataBase/FOM_Captures/{simulation_name}.png"
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

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def UpdateProjectParametersFOM(parameters, mu=None, mesh_file_name=None):
    angle_of_attack        = mu[0]
    mach_infinity          = mu[1]
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(angle_of_attack)
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach_infinity)
    parameters["output_processes"]["vtk_output"][0]["Parameters"]["output_path"].SetString(f'DataBase/FOM_Results/{angle_of_attack}, {mach_infinity}')
    parameters["modelers"][0]["parameters"]["input_filename"].SetString(f'{mesh_file_name}.med')
    return parameters

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def _GetAnalysisStageClass(parameters):
    analysis_stage_module_name = parameters["analysis_stage"].GetString()
    analysis_stage_class_name = analysis_stage_module_name.split('.')[-1]
    analysis_stage_class_name = ''.join(x.title() for x in analysis_stage_class_name.split('_'))
    analysis_stage_module = importlib.import_module(analysis_stage_module_name)
    analysis_stage_class = getattr(analysis_stage_module, analysis_stage_class_name)
    return analysis_stage_class

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
    if os.path.exists("mu_train.npy") and os.path.exists("mu_test.npy"):
        archivo = 'mu_train.npy'
        mu_train = np.load(archivo)
        archivo = 'mu_test.npy'
        mu_test = np.load(archivo)
    elif os.path.exists("mu_train.npy"):
        archivo = 'mu_train.npy'
        mu_train = np.load(archivo)
        mu_test = []
    elif os.path.exists("mu_test.npy"):
        archivo = 'mu_test.npy'
        mu_test = np.load(archivo)
        mu_train = []
    return np.array(mu_train), np.array(mu_test)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# plot parameters
#
def plot_mu_values(mu_train, mu_test):
    if len(mu_train) > 0: plt.plot(mu_train[:,1], mu_train[:,0], 'bs', label="Train Values")
    if len(mu_test) > 0: plt.plot(mu_test[:,1], mu_test[:,0], 'ro', label="Test Values")
    plt.title('Mu Values')
    plt.ylabel('Alpha')
    plt.xlabel('Mach')
    plt.grid(True)
    plt.legend(bbox_to_anchor=(.85, 1.03, 1., .102), loc='upper left', borderaxespad=0.)
    plt.savefig("MuValues.png")
    plt.close('all')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# PEBL clustering
# 
def error_estimation(mu_train, mu_test):
    if len(mu_train) > 0:
        approximation_error = 0.0
        FOM_model = []; RBF_model = []
        for mu in mu_train:
            FOM_model.append(np.load(f'DataBase/Snapshots/{mu[0]}, {mu[1]}.npy'))
            RBF_model.append(np.load(f"model_mesh_0/Interpolated_data/{mu[0]}, {mu[1]}.npy").T)
        FOM_model = np.block(FOM_model)
        RBF_model = np.block(RBF_model)
        training_approximation_error = np.linalg.norm(FOM_model - RBF_model)/np.linalg.norm(FOM_model)
        print(f'PEBL CLUSTERING: RBF training approximation error: {training_approximation_error}')

    if len(mu_test) > 0:
        approximation_error = 0.0
        FOM_model = []; RBF_model_interpolation = []
        for mu in mu_test:
            FOM_model.append(np.load(f'DataBase/Snapshots/{mu[0]}, {mu[1]}.npy'))
            RBF_model_interpolation.append(np.load(f"model_mesh_0/Interpolated_data/{mu[0]}, {mu[1]}.npy").T)
        FOM_model = np.block(FOM_model)
        RBF_model_interpolation = np.block(RBF_model_interpolation)
        approximation_error = np.linalg.norm(FOM_model - RBF_model_interpolation)/np.linalg.norm(FOM_model)
        print(f'PEBL CLUSTERING: RBF interpolation approximation error: {approximation_error}')

def E_p(u, c):
    """
    c: direction vector onto which to project.
    u: vector or collection of column vectors to project onto the direction of c.
    """
    c = c.reshape(-1, 1)
    if len(u.shape) == 1:
        u = u.reshape(-1, 1)
    projection_coefficients = (u.T @ c) / (c.T @ c)
    projection_of_u_onto_c = projection_coefficients.T * c
    projection_error = np.linalg.norm(u - projection_of_u_onto_c, axis=0) / np.linalg.norm(u, axis=0)
    return projection_error

def ObtainBasis(Snapshots, truncation_tolerance=0):
        U,_,_= truncated_svd(Snapshots,truncation_tolerance)
        return U

def truncated_svd(Matrix, epsilon=0):
    M,N=np.shape(Matrix)
    dimMATRIX = max(M,N)
    U, s, V = np.linalg.svd(Matrix, full_matrices=False) #U --> M xN, V --> N x N
    V = V.T
    tol = dimMATRIX*np.finfo(float).eps*max(s)/2
    R = np.sum(s > tol)  # Definition of numerical rank
    if epsilon == 0:
        K = R
    else:
        SingVsq = np.multiply(s,s)
        SingVsq.sort()
        normEf2 = np.sqrt(np.cumsum(SingVsq))
        epsilon = epsilon*normEf2[-1] #relative tolerance
        T = (sum(normEf2<epsilon))
        K = len(s)-T
    K = min(R,K)
    return U[:, :K], s[:K], V[:, :K]

class PEBL:
    def __init__(self, bisection_tolerance=0.15,  POD_tolerance=1e-3):
        self.bisection_tolerance = bisection_tolerance
        self.POD_tolerance = POD_tolerance

    def fit(self, Snapshots):
        self.Snapshots = Snapshots
        self.Stage_1()
        self.Stage_2()
        return self

    def Stage_1(self):
        #stage 1, generation of bisection tree with accuracy 'bisection_tolerance'
        max_index = np.argmax( np.linalg.norm(self.Snapshots, axis=0) )
        first_snapshot = self.Snapshots[:,max_index]
        self.Tree = Node([first_snapshot,np.arange(0,self.Snapshots.shape[1], 1, dtype=int)])
        bisect_flag = True
        while bisect_flag == True:
            bisect_flag = False
            for leaf in self.Tree.leaves():
                errors = E_p(self.Snapshots[:,leaf.val[1]], leaf.val[0])
                max_error = max(errors)
                print(f'Clustering local max error {max_error}')
                if max_error > self.bisection_tolerance:
                    bisect_flag = True
                    #find next anchor point
                    max_index = np.argmax(errors)
                    c_new = self.Snapshots[:,leaf.val[1]][:,max_index]
                    new_errors = E_p(self.Snapshots[:,leaf.val[1]], c_new)
                    indexes_left = np.where( errors <= new_errors)
                    indexes_right = np.where( errors > new_errors)
                    #divide the snapshots among the two children
                    leaf.left =  Node([leaf.val[0], leaf.val[1][indexes_left[0]]])
                    leaf.right = Node([c_new, leaf.val[1][indexes_right[0]]])
                    leaf.val[1] = None

    def Stage_2(self):
        #stage 2, generation of local POD bases'
        index = 0
        for leaf in self.Tree.leaves():
            Phi_i = ObtainBasis(self.Snapshots[:,leaf.val[1]], self.POD_tolerance)
            leaf.val.append(Phi_i)
            leaf.val.append(index)
            index+=1

    def predict(self, u):
        current_node = self.Tree
        while not current_node.is_leaf():
            if E_p(u, current_node.left.val[0]) < E_p(u, current_node.right.val[0]):
                current_node = current_node.left
            else:
                current_node = current_node.right
        return current_node.val[3]

class Node:
    def __init__(self, val):
        self.left = None
        self.right = None
        self.val = val

    def leaves(self):
        current_nodes = [self]
        leaves = []

        while len(current_nodes) > 0:
            next_nodes = []
            for node in current_nodes:
                if node.left is None and node.right is None:
                    leaves.append(node)
                    continue
                if node.left is not None:
                    next_nodes.append(node.left)
                if node.right is not None:
                    next_nodes.append(node.right)
            current_nodes = next_nodes
        return leaves

    def is_leaf(self):
        return self.left is None and self.right is None
    
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# get multiple parameters by Halton or LatinHypercube methods
#
def get_multiple_parameters(number_train_values=0, number_test_values=0, angle=[], mach=[], method='Halton'):
    from scipy.stats import qmc                 
    if method == 'Halton':
        sampler = qmc.Halton(d=2)
    elif method == 'LatinHypercube':
        sampler = qmc.LatinHypercube(d=2)
    mu_train = []; mu_test = []
    if number_train_values > 0:
        sample = sampler.random(number_train_values)
        values = qmc.scale(sample, [angle[0],mach[0]], [angle[1],mach[1]])
        for i in range(number_train_values):
            #Angle of attack , Mach infinit
            mu_train.append([values[i,0], values[i,1]])
    if number_test_values > 0:
        sample = sampler.random(number_test_values)
        values = qmc.scale(sample, [angle[0],mach[0]], [angle[1],mach[1]])
        for i in range(number_test_values):
            #Angle of attack , Mach infinit
            mu_test.append([values[i,0], values[i,1]])
    return mu_train, mu_test

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def GetRomManagerParameters():
    general_rom_manager_parameters = KratosMultiphysics.Parameters("""{
            "rom_stages_to_train" : ["ROM"],      // ["ROM","HROM"]
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

    update_parameters        = True
    update_bases             = False
    update_clusters          = False
    plot_custering_info      = True
    plot_clusters            = True
    plot_clusters_prediction = True
    ###############################
    # MESH LIST
    meshes_list = ['model_mesh_0']
    ###############################
    # PARAMETERS SETTINGS
    number_of_mu_train = 5
    number_of_mu_test  = 1
    mach_range  = [ 0.70, 0.75]
    angle_range = [ 1.00, 2.00]
    ###############################
    # CLUSTERING SETTINGS
    bisection_tolerance = 0.15
    POD_tolerance       = 1e-12
    ###############################
    
    if update_parameters:
        mu_train, mu_test = get_multiple_parameters(
            number_train_values = number_of_mu_train,
            number_test_values  = number_of_mu_test , 
            angle               = angle_range       , 
            mach                = mach_range        , 
            method              = 'Halton'          )
        
        for mesh_file in meshes_list:
            if not os.path.exists(f"DataBase"):
                os.mkdir(f"DataBase")
            else:
                KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting(f"DataBase")
                os.mkdir(f"DataBase")
            if not os.path.exists(f"DataBase/Snapshots"):
                os.mkdir(f"DataBase/Snapshots")
            if not os.path.exists(f"DataBase/Snapshots_not_converged_steps"):
                os.mkdir(f"DataBase/Snapshots_not_converged_steps")
            mu_list = []
            with os.scandir(f'DataBase/Snapshots') as files:
                files = [file.name for file in files if file.is_file() and file.name.endswith('.npy')]
            print(f'Initial number of snapshots: {len(files)}')
            mu_aux_list = list(mu_train) + list(mu_test)
            for mu in mu_aux_list:
                file = f'{mu[0]}, {mu[1]}'
                if (   not os.path.exists(f"DataBase/FOM_Captures/{file}.png")
                    or not os.path.exists(f"DataBase/FOM_Results/{file}")
                    or not os.path.exists(f"DataBase/FOM_Skin_Data/{file}.dat")
                    or not os.path.exists(f"DataBase/Snapshots/{file}.npy")
                    or not os.path.exists(f"DataBase/Snapshots_not_converged_steps/{file}.npy")):
                    mu_list.append([np.double(mu[0]), np.double(mu[1])])

            LaunchFOMSimulations(mu_list, mesh_file)

            snapshotsMatrix = []
            aux_snapshotsMatrix = []
            mu_aux_list = list(mu_train)
            for mu in mu_aux_list:
                file = f'{mu[0]}, {mu[1]}.npy'
                snapshotsMatrix.append(np.load(f'DataBase/Snapshots/{file}'))
                aux_snapshotsMatrix.append(np.load(f'DataBase/Snapshots_not_converged_steps/{file}'))
            snapshotsMatrix = np.block(snapshotsMatrix)
            aux_snapshotsMatrix = np.block(aux_snapshotsMatrix)
            print(snapshotsMatrix.shape)
            print(aux_snapshotsMatrix.shape)
            np.save(f'SnapshotsMatrix_conv', snapshotsMatrix)
            np.save(f'SnapshotsMatrix_not_conv', aux_snapshotsMatrix)
            KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting(f"rom_data")
            KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting(f"__pycache__")

    else:
            mu_train, mu_test = load_mu_parameters()
    
    plot_mu_values(np.array(mu_train), np.array(mu_test))
    save_mu_parameters(mu_train, 'mu_train', mu_test, 'mu_test')

    if update_parameters: update_clusters = True
    if update_clusters:
        #### CLUSTERING DATA
        #################################################################################################
        parameters = np.array(mu_train)
        snapshots = np.load(f'SnapshotsMatrix_conv.npy')

        if plot_custering_info: print(f'PEBL CLUSTERING: Snapshots Matrix shape: {snapshots.shape}')
        if plot_custering_info: print(f'PEBL CLUSTERING: Parameters shape: {parameters.shape}')

        #### RBF TRAINING
        #################################################################################################
        db = Database(parameters, snapshots.T)
        pod = POD()
        rbf = RBF()
        rom = ROM(db, pod, rbf).fit()
        
        #### PEBL CLUSTERING 
        #################################################################################################
        pebl_object = PEBL(bisection_tolerance=bisection_tolerance,  POD_tolerance=POD_tolerance).fit(snapshots)
        leaves = pebl_object.Tree.leaves()
        mu_train_aux = np.array(mu_train)
        mu_train_with_leaves = []
        leaf_id = 0
        for leaf in leaves:
            for i in range(len(leaf.val[1])):
                mu_train_with_leaves.append([mu_train_aux[leaf.val[1]][i,0], mu_train_aux[leaf.val[1]][i,1], leaf_id])
            leaf_id += 1
        if not os.path.exists(f"model_mesh_0"):
            os.mkdir(f"model_mesh_0")
        np.save(f'model_mesh_0/mu_train_with_indexes', mu_train_with_leaves)

        if plot_custering_info: print(f'PEBL CLUSTERING: Number of leaves (clusters): {len(leaves)}')

        training_solutions_list = [rom.predict([element]).snapshots_matrix for element in mu_train]
        if not os.path.exists(f"model_mesh_0/Interpolated_data"):
            os.mkdir(f"model_mesh_0/Interpolated_data")
        for i, solution in enumerate(training_solutions_list):
            np.save(f"model_mesh_0/Interpolated_data/{mu_train[i][0]}, {mu_train[i][1]}.npy", solution)
        
        if plot_clusters:
            #### PLOT LEAVES / CLUSTERS
            #################################################################################################
            fig = plt.figure()
            fig.set_figwidth(24.0)
            fig.set_figheight(16.0)
            leaf_id = 0
            for leaf in leaves:
                if plot_custering_info: print("PEBL CLUSTERING: Number of cases in leaf "+str(leaf_id)+": ", len(leaf.val[1]))
                fig = plt.scatter(mu_train_aux[leaf.val[1]][:,1], mu_train_aux[leaf.val[1]][:,0], label="Cluster "+str(leaf_id))
                fig = plt.scatter(np.mean(mu_train_aux[leaf.val[1]][:,1]), np.mean(mu_train_aux[leaf.val[1]][:,0]), c='k',marker=f'${leaf_id}$', s= 150)
                leaf_id += 1
            fig = plt.title("PEBL clustering")
            fig = plt.ylabel('Alpha')
            fig = plt.xlabel('Mach')
            fig = plt.legend(bbox_to_anchor=(1.004, 0.9, 1.0, 0.102), loc='upper left', borderaxespad=0.)
            fig = plt.savefig(f'model_mesh_0/PEBL clustering.png',bbox_inches='tight' )
            fig = plt.close('all')

        if len(mu_test) > 0:
            #### SOLUTION AND CLUSTER PREDICTION TEST
            #################################################################################################
            if not os.path.exists(f"model_mesh_0/Interpolated_data"):
                os.mkdir(f"model_mesh_0/Interpolated_data")
            interpolated_solutions_list = [rom.predict([element]).snapshots_matrix for element in mu_test]
            predicted_indexes_list = [pebl_object.predict(solution.T) for solution in interpolated_solutions_list]
            for i, solution in enumerate(interpolated_solutions_list):
                np.save(f"model_mesh_0/Interpolated_data/{mu_test[i][0]}, {mu_test[i][1]}.npy", solution)
            np.save(f'model_mesh_0/predicted_indexes', predicted_indexes_list)

            if plot_clusters_prediction:
                #### PLOT LEAVES / CLUSTERS PREDICTION
                #################################################################################################
                if not os.path.exists(f"model_mesh_0/Clustering_mu_test_plots"):
                    os.mkdir(f"model_mesh_0/Clustering_mu_test_plots")
                for i, predicted_index in enumerate(predicted_indexes_list):
                    if plot_custering_info: print(f'Test {i}: {mu_test[i][1]} {mu_test[i][0]}, predicted index: {predicted_index}')
                    #### PLOT CLUSTERS PREDICTION
                    fig = plt.figure()
                    fig.set_figwidth(24.0)
                    fig.set_figheight(16.0)
                    leaf_id = 0
                    for leaf in leaves:
                        if leaf_id == predicted_index:
                            fig = plt.scatter(mu_train_aux[leaf.val[1]][:,1], mu_train_aux[leaf.val[1]][:,0], label="Cluster "+str(leaf_id))
                            fig = plt.scatter(np.mean(mu_train_aux[leaf.val[1]][:,1]), np.mean(mu_train_aux[leaf.val[1]][:,0]), c='k',marker=f'${leaf_id}$', s= 150)
                        else:
                            fig = plt.scatter(mu_train_aux[leaf.val[1]][:,1], mu_train_aux[leaf.val[1]][:,0], label="Cluster "+str(leaf_id), alpha=0.2)
                            fig = plt.scatter(np.mean(mu_train_aux[leaf.val[1]][:,1]), np.mean(mu_train_aux[leaf.val[1]][:,0]), c='k',marker=f'${leaf_id}$', alpha=0.2, s= 150)
                        leaf_id += 1
                    fig = plt.scatter(mu_test[i][1], mu_test[i][0], c='k',marker="s", s=150)
                    fig = plt.text(mu_test[i][1], mu_test[i][0]-0.05, 
                                    f'Mu test {i}', ha='center',  
                                    bbox=dict(boxstyle="round", ec=(1., 0.5, 0.5),fc=(1., 0.8, 0.8),))
                    fig = plt.title("PEBL online prediction")
                    fig = plt.ylabel('Alpha')
                    fig = plt.xlabel('Mach')
                    fig = plt.legend(bbox_to_anchor=(1.004, 0.9, 1.0, 0.102), loc='upper left', borderaxespad=0.)
                    fig = plt.savefig(f'model_mesh_0/Clustering_mu_test_plots/PEBL prediction {i} {mu_test[i][1]} {mu_test[i][0]}.png',bbox_inches='tight' )
                    fig = plt.close('all')
        
        if plot_custering_info: error_estimation(mu_train, mu_test)   

    if update_parameters or update_clusters: update_bases = True
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    general_rom_manager_parameters = GetRomManagerParameters()
    project_parameters_name = "ROMProjectParameters.json"
    rom_manager = LocalRomManager(project_parameters_name,general_rom_manager_parameters,CustomizeSimulation,UpdateProjectParameters)
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    mesh_errors = []
    for mesh in meshes_list:

        if not os.path.exists(f"{mesh}"):
            os.mkdir(f"{mesh}")

        if os.path.exists(f'{mesh}/mu_train_with_indexes.npy'):  
            start_time = time.time()

            mu_train_with_indexes = np.load(f'{mesh}/mu_train_with_indexes.npy')
            mu_train_errors = rom_manager.Fit(mu_train, mu_train_with_indexes, mesh, update_bases)

            train_exe_time = time.time() - start_time
            print('Executing ROM Manager took ' + str(round(train_exe_time, 2)) + ' sec')

            _,train_error,_,_ = rom_manager.PrintErrors()

            mesh_errors.append([mesh, 'Train', round(train_exe_time, 2), train_error])
            print(f'Mesh name | Execution time | Approximation error FOM vs ROM')
            print(mesh_errors)
        else: 
            mu_train_errors = []

        if os.path.exists(f'{mesh}/predicted_indexes.npy'): 
            start_time = time.time()
            
            predicted_indexes     = np.load(f'{mesh}/predicted_indexes.npy')
            mu_test_errors  = rom_manager.Test(mu_test, mesh, predicted_indexes)

            test_exe_time = time.time() - start_time
            print('Executing ROM Manager took ' + str(round(test_exe_time, 2)) + ' sec')

            _,_,_,test_error = rom_manager.PrintErrors()

            mesh_errors.append([mesh, 'Test', round(test_exe_time, 2), test_error])
            print(f'    Mesh name    |    Execution time    |    Approximation error FOM vs ROM')
            print(mesh_errors[0])
            print(mesh_errors[1])
        else:
            mu_test_errors = []

        save_mu_parameters(mu_train_errors, f'{mesh}/mu_train_errors', mu_test_errors, f'{mesh}/mu_test_errors')