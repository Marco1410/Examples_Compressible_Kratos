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

def load_mu_parameters():
    if os.path.exists("../mu_parameters/mu_train.npy") and os.path.exists("../mu_parameters/mu_test.npy"):
        archivo = '../mu_parameters/mu_train.npy'
        mu_train = np.load(archivo)
        archivo = '../mu_parameters/mu_test.npy'
        mu_test = np.load(archivo)
    elif os.path.exists("../mu_parameters/mu_train.npy"):
        archivo = '../mu_parameters/mu_train.npy'
        mu_train = np.load(archivo)
        mu_test = []
    elif os.path.exists("../mu_parameters/mu_test.npy"):
        archivo = '../mu_parameters/mu_test.npy'
        mu_test = np.load(archivo)
        mu_train = []
    return np.array(mu_train), np.array(mu_test)

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
    with open('../ProjectParametersFiles/FomProjectParameters.json','r') as parameter_file:
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

def UpdateProjectParameters(parameters, mu=None, mesh_file_name=None):
    angle_of_attack        = mu[0]
    mach_infinity          = mu[1]
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(angle_of_attack)
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach_infinity)
    parameters["output_processes"]["vtk_output"][0]["Parameters"]["output_path"].SetString(f'{mesh_file_name}/FOM_Results/{angle_of_attack}, {mach_infinity}')
    parameters["modelers"][0]["parameters"]["input_filename"].SetString(f'../salome_mesh_files/{mesh_file_name}.med')
    return parameters

def _GetAnalysisStageClass(parameters):
    analysis_stage_module_name = parameters["analysis_stage"].GetString()
    analysis_stage_class_name = analysis_stage_module_name.split('.')[-1]
    analysis_stage_class_name = ''.join(x.title() for x in analysis_stage_class_name.split('_'))
    analysis_stage_module = importlib.import_module(analysis_stage_module_name)
    analysis_stage_class = getattr(analysis_stage_module, analysis_stage_class_name)
    return analysis_stage_class

if __name__ == '__main__':

    start_time = time.time()

    with os.scandir('../salome_mesh_files') as files:
        mesh_files = [file.name.removesuffix(".med") for file in files if file.is_file() and file.name.endswith('.med') and 'SMESH' not in file.name]
    mesh_files.sort()

    mu_train, mu_test = load_mu_parameters()

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