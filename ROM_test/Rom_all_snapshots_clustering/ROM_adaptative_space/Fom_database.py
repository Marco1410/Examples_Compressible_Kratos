import os
import openpyxl
import importlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import time as time
from parameters import * # type: ignore
import KratosMultiphysics
import KratosMultiphysics.kratos_utilities
from KratosMultiphysics.RomApplication.calculate_rom_basis_output_process import CalculateRomBasisOutputProcess

def LaunchFOMSimulations(mu_test):

    with open('ProjectParameters.json','r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    parameters["output_processes"].AddEmptyArray("rom_output")
    parameters["output_processes"]["rom_output"].Append(_GetDefaulRomBasisOutputParameters())

    for Id, mu in enumerate(mu_test):
        fom_info_steps_list = []
        model = KratosMultiphysics.Model()
        parameters_copy = UpdateProjectParameters(parameters.Clone(), mu)
        analysis_stage_class = _GetAnalysisStageClass(parameters_copy)
        simulation = CustomizeSimulation(analysis_stage_class, model, parameters_copy)

        start_time = time.time()
        simulation.Run()
        exe_time = time.time() - start_time

        fom_info_steps_list.append([Id,
                                    mu[0],
                                    mu[1], 
                                    simulation.model["MainModelPart"].ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER],
                                    simulation.model["MainModelPart"].ProcessInfo[KratosMultiphysics.RESIDUAL_NORM],
                                    round(exe_time, 2)])
        
        for process in simulation._GetListOfOutputProcesses():
            if isinstance(process, CalculateRomBasisOutputProcess):
                BasisOutputProcess = process

        np.save(f'DataBase/Snapshots/{mu[0]}, {mu[1]}',BasisOutputProcess._GetSnapshotsMatrix()) 
        np.save(f'DataBase/Snapshots_not_converged/{mu[0]}, {mu[1]}', np.array(simulation._GetSolver()._GetSolutionStrategy().GetIntermediateSolutionsMatrix()))

        if os.path.exists(f"fom_data.xlsx"):
            wb = openpyxl.load_workbook(f"fom_data.xlsx")
            hoja = wb.active
            for item in fom_info_steps_list:
                hoja.append(item)
            wb.save(f'fom_data.xlsx')
        else:
            wb = openpyxl.Workbook()
            hoja = wb.active
            hoja.append(('Id', 'Angle', 'Mach', 'NL iterations', 'Residual norm', 'Time [sec]'))
            for item in fom_info_steps_list:
                hoja.append(item)
            wb.save(f'fom_data.xlsx')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('rom_data')

def CustomizeSimulation(cls, global_model, parameters):

    class CustomSimulation(cls):
        def __init__(self, model,project_parameters):
            super().__init__(model,project_parameters)

        def Initialize(self):
            super().Initialize()

        def InitializeSolutionStep(self):
            super().InitializeSolutionStep()

        def FinalizeSolutionStep(self):
            super().FinalizeSolutionStep()
            nametype = parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].GetString()
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
            fig = plt.plot(x, cp, 'xr', markersize = 2.0, label = 'FOM')
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

def UpdateProjectParameters(parameters, mu=None):
    angle_of_attack        = mu[0]
    mach_infinity          = mu[1]
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(angle_of_attack)
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach_infinity)
    parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString(f'DataBase/FOM_Results/{angle_of_attack}, {mach_infinity}')
    return parameters

def _GetAnalysisStageClass(parameters):
    analysis_stage_module_name = parameters["analysis_stage"].GetString()
    analysis_stage_class_name = analysis_stage_module_name.split('.')[-1]
    analysis_stage_class_name = ''.join(x.title() for x in analysis_stage_class_name.split('_'))
    analysis_stage_module = importlib.import_module(analysis_stage_module_name)
    analysis_stage_class = getattr(analysis_stage_module, analysis_stage_class_name)
    return analysis_stage_class

def _GetDefaulRomBasisOutputParameters():
    return KratosMultiphysics.Parameters("""{
            "python_module" : "calculate_rom_basis_output_process",
            "kratos_module" : "KratosMultiphysics.RomApplication",
            "process_name"  : "CalculateRomBasisOutputProcess",
            "help"          : "This process should write the Rom basis",
            "Parameters"    :
            {
                "model_part_name": "MainModelPart",
                "rom_manager" : false,      // set to false for manual manipulation of ROM via flags in the RomParameters
                "snapshots_control_type": "step",
                "snapshots_interval": 1.0,
                "nodal_unknowns":  ["VELOCITY_POTENTIAL","AUXILIARY_VELOCITY_POTENTIAL"],
                "rom_basis_output_format": "numpy",
                "rom_basis_output_name": "rom_bases",
                "rom_basis_output_folder": "rom_data",
                "svd_truncation_tolerance": 0
            }
        }""")
    
if __name__ == '__main__':

    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('DataBase')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('fom_data.xlsx')
    os.mkdir(f"DataBase/Snapshots")
    os.mkdir(f"DataBase/Snapshots_not_converged")
    os.mkdir(f"DataBase/FOM_Results")
    os.mkdir(f"DataBase/FOM_Skin_Data")
    os.mkdir(f"DataBase/FOM_Captures")

    mu_train,_, mu_test,_ = load_mu_parameters() # type: ignore

    start_time = time.time()
    LaunchFOMSimulations(list(mu_train) + list(mu_test))
    exe_time = time.time() - start_time

    print(f' Executing took {round(exe_time, 2)} sec')

    snapshotsMatrix = []
    snapshotsMatrix_nc = []
    for mu in list(mu_train):
        file = f'{mu[0]}, {mu[1]}.npy'
        snapshotsMatrix.append(np.load(f'DataBase/Snapshots/{file}'))
        snapshotsMatrix_nc.append(np.load(f'DataBase/Snapshots_not_converged/{file}'))
    snapshotsMatrix = np.block(snapshotsMatrix)
    snapshotsMatrix_nc = np.block(snapshotsMatrix_nc)
    print(f'SnapshotsMatrix shape: {snapshotsMatrix.shape}')
    print(f'SnapshotsMatrix not converged step shape: {snapshotsMatrix_nc.shape}')
    np.save(f'DataBase/SnapshotsMatrix_conv', snapshotsMatrix)
    np.save(f'DataBase/SnapshotsMatrix_not_conv', snapshotsMatrix_nc)
