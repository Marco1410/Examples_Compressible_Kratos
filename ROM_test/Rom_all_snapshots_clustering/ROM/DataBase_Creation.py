import os
import openpyxl
import importlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import time as time
import KratosMultiphysics
from KratosMultiphysics.RomApplication.calculate_rom_basis_output_process import CalculateRomBasisOutputProcess

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
#
# save parameters
#
def save_mu_parameters(mu_train, name1, mu_test, name2):
    if len(mu_train) > 0:
        np.save(f'{name1}',mu_train)
    if len(mu_test) > 0:
        np.save(f'{name2}',mu_test)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# save parameters
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
# save parameters
#
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
# save parameters
#
def LaunchSimulations(mu_test, mesh_file_name):
    if not os.path.exists(f"{mesh_file_name}/Snapshots"):
        os.mkdir(f"{mesh_file_name}/Snapshots")
    if not os.path.exists(f"{mesh_file_name}/Snapshots_not_converged_steps"):
        os.mkdir(f"{mesh_file_name}/Snapshots_not_converged_steps")
    if not os.path.exists(f"{mesh_file_name}/FOM_Results"):
        os.mkdir(f"{mesh_file_name}/FOM_Results")
    if not os.path.exists(f"{mesh_file_name}/FOM_Skin_Data"):
        os.mkdir(f"{mesh_file_name}/FOM_Skin_Data")
    if not os.path.exists(f"{mesh_file_name}/FOM_Captures"):
        os.mkdir(f"{mesh_file_name}/FOM_Captures")
    with open('FomProjectParameters.json','r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    for Id, mu in enumerate(mu_test):
        fom_info_steps_list = []
        model = KratosMultiphysics.Model()
        parameters_copy = UpdateProjectParameters(parameters.Clone(), mu, mesh_file_name)
        analysis_stage_class = _GetAnalysisStageClass(parameters_copy)
        simulation = CustomizeSimulation(analysis_stage_class, model, parameters_copy, mesh_file_name)
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
        np.save(f'{mesh_file_name}/Snapshots/{mu[0]}, {mu[1]}',BasisOutputProcess._GetSnapshotsMatrix()) 
        np.save(f'{mesh_file_name}/Snapshots_not_converged_steps/{mu[0]}, {mu[1]}', np.array(simulation._GetSolver()._GetSolutionStrategy().GetIntermediateSolutionsMatrix()))

        if os.path.exists(f"{mesh_file_name}/fom_data.xlsx"):
            wb = openpyxl.load_workbook(f"{mesh_file_name}/fom_data.xlsx")
            hoja = wb.active
            for item in fom_info_steps_list:
                hoja.append(item)
            wb.save(f'{mesh_file_name}/fom_data.xlsx')
        else:
            wb = openpyxl.Workbook()
            hoja = wb.active
            hoja.append(('Id', 'Angle', 'Mach', 'NL iterations', 'Residual norm'))
            for item in fom_info_steps_list:
                hoja.append(item)
            wb.save(f'{mesh_file_name}/fom_data.xlsx')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# save parameters
#
def CustomizeSimulation(cls, global_model, parameters, mesh_file_name):

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
            simulation_name = nametype.removeprefix(f"{mesh_file_name}/FOM_Results/")
            skin_data_filename = f"{mesh_file_name}/FOM_Skin_Data/{simulation_name}.dat"
            capture_filename = f"{mesh_file_name}/FOM_Captures/{simulation_name}.png"
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
#
# save parameters
#
def UpdateProjectParameters(parameters, mu=None, mesh_file_name=None):
    angle_of_attack        = mu[0]
    mach_infinity          = mu[1]
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(angle_of_attack)
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach_infinity)
    parameters["output_processes"]["vtk_output"][0]["Parameters"]["output_path"].SetString(f'{mesh_file_name}/FOM_Results/{angle_of_attack}, {mach_infinity}')
    parameters["modelers"][0]["parameters"]["input_filename"].SetString(f'{mesh_file_name}.med')
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

if __name__ == '__main__':

    start_time = time.time()

    mesh_files = ['model_mesh_0']

    NumberofMuTrain = 1000
    NumberOfMuTest  = 1000

    # Definir rango de valores de mach y angulo de ataque
    mach_range  = [ 0.65, 0.75]
    angle_range = [ 0.00, 3.00]

    mu_train, mu_test = get_multiple_parameters(number_train_values = NumberofMuTrain,
                        number_test_values  = NumberOfMuTest, 
                        angle  = angle_range, 
                        mach   = mach_range, 
                        method = 'Halton')
    
    # mu_train, mu_test = load_mu_parameters()
    
    mu_train = np.array(mu_train)
    mu_test  = np.array(mu_test)
    
    plot_mu_values(mu_train, mu_test)

    save_mu_parameters(mu_train, 'mu_train', mu_test, 'mu_test')

    for mesh_file in mesh_files:
        if not os.path.exists(f"{mesh_file}"):
            os.mkdir(f"{mesh_file}")
        if not os.path.exists(f"{mesh_file}/Snapshots"):
            os.mkdir(f"{mesh_file}/Snapshots")
        if not os.path.exists(f"{mesh_file}/Snapshots_not_converged_steps"):
            os.mkdir(f"{mesh_file}/Snapshots_not_converged_steps")
        mu_list = []
        with os.scandir(f'{mesh_file}/Snapshots') as files:
            files = [file.name for file in files if file.is_file() and file.name.endswith('.npy')]
        print(f'Initial number of snapshots: {len(files)}')
        mu_aux_list = list(mu_train) + list(mu_test)
        for mu in mu_aux_list:
            file = f'{mu[0]}, {mu[1]}'
            if (   not os.path.exists(f"{mesh_file}/FOM_Captures/{file}.png")
                or not os.path.exists(f"{mesh_file}/FOM_Results/{file}")
                or not os.path.exists(f"{mesh_file}/FOM_Skin_Data/{file}.dat")
                or not os.path.exists(f"{mesh_file}/Snapshots/{file}.npy")
                or not os.path.exists(f"{mesh_file}/Snapshots_not_converged_steps/{file}.npy")):
                mu_list.append([np.double(mu[0]), np.double(mu[1])])

        LaunchSimulations(mu_list, mesh_file)

        with os.scandir(f'{mesh_file}/Snapshots') as files:
            files = [file.name for file in files if file.is_file() and file.name.endswith('.npy')]
        print(f'Final number of snapshots: {len(files)}')

        snapshotsMatrix = []
        aux_snapshotsMatrix = []
        mu_aux_list = list(mu_train)
        for mu in mu_aux_list:
            file = f'{mu[0]}, {mu[1]}.npy'
            snapshotsMatrix.append(np.load(f'{mesh_file}/Snapshots/{file}'))
            aux_snapshotsMatrix.append(np.load(f'{mesh_file}/Snapshots_not_converged_steps/{file}'))
        snapshotsMatrix = np.block(snapshotsMatrix)
        aux_snapshotsMatrix = np.block(aux_snapshotsMatrix)
        print(snapshotsMatrix.shape)
        print(aux_snapshotsMatrix.shape)
        np.save(f'{mesh_file}/SnapshotsMatrix_conv', snapshotsMatrix)
        np.save(f'{mesh_file}/SnapshotsMatrix_not_conv', aux_snapshotsMatrix)

    exe_time = time.time() - start_time

    print(f' Executing took {round(exe_time, 2)} sec')