import os
import importlib
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
import KratosMultiphysics
import KratosMultiphysics.kratos_utilities
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
from KratosMultiphysics.gid_output_process import GiDOutputProcess


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# save parameters
#
def LaunchFakeSimulation(data_set, mu):
    with open('ProjectParameters.json','r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
        model = KratosMultiphysics.Model()
        parameters_copy = UpdateProjectParameters(parameters.Clone(), mu)
        analysis_stage_class = _GetAnalysisStageClass(parameters_copy)
        simulation = CustomizeSimulation(analysis_stage_class, model, parameters_copy, data_set)
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

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# save parameters
#
def CustomizeSimulation(cls, global_model, parameters, data_set):

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
            nametype = parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].GetString()
            simulation_name = nametype.removeprefix(f"RBF_Results/")
            skin_data_filename = f"RBF_Skin_Data/{simulation_name}.dat"
            fout = open(skin_data_filename,'w')
            modelpart = self.model["MainModelPart.Body2D_Body"]
            for node in modelpart.Nodes:
                x = node.X ; y = node.Y ; z = node.Z
                cp = node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
                fout.write("%s %s %s %s\n" %(x,y,z,cp))
            fout.close()

    return CustomSimulation(global_model, parameters)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# save parameters
#
def UpdateProjectParameters(parameters, mu=None):
    angle_of_attack        = mu[0]
    mach_infinity          = mu[1]
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(angle_of_attack)
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach_infinity)
    parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString(f'RBF_Results/{angle_of_attack}, {mach_infinity}')
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

    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Train_fom_rom_rbf_errors') 
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Test_fom_rom_rbf_errors') 
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('RBF_Skin_Data') 
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('RBF_Results') 

    train = True
    test  = True

    if os.path.exists('mu_train.npy') and train:
        mu_train = np.load('mu_train.npy')

        for mu in mu_train:

            if not os.path.exists('Train_fom_rom_rbf_errors'):
                os.mkdir('Train_fom_rom_rbf_errors')

            if not os.path.exists('RBF_Skin_Data'):
                os.mkdir('RBF_Skin_Data')

            capture_name = f"Train_fom_rom_rbf_errors/Fit_{mu[0]}, {mu[1]}.png"

            #### ONLINE CP PLOT
            ######################################################################
            cp_min = 0
            cp_max = 0
            fig = plt.figure()
            fig.set_figwidth(12.0)
            fig.set_figheight(8.0)

            #### FOM ####
            fom_skin_data_filename = f"DataBase/FOM_Skin_Data/{mu[0]}, {mu[1]}.dat"
            if os.path.exists(fom_skin_data_filename):
                x_fom  = np.loadtxt(fom_skin_data_filename, usecols=(0,))
                cp_fom = np.loadtxt(fom_skin_data_filename, usecols=(3,))
                fig = plt.plot(x_fom, cp_fom, 'ob', markersize = 2.0, label = 'FOM')

            #### ROM ####
            rom_skin_data_filename = f"Skin_Data/ROM_Fit{mu[0]}, {mu[1]}.dat"
            if os.path.exists(rom_skin_data_filename):
                fom = np.load(f'DataBase/Snapshots/{mu[0]}, {mu[1]}.npy')
                rom = np.load(f'RomDataBase/{mu[0]}, {mu[1]}.npy')
                x_rom  = np.loadtxt(rom_skin_data_filename, usecols=(0,))
                cp_rom = np.loadtxt(rom_skin_data_filename, usecols=(3,))
                fig = plt.plot(x_rom, cp_rom, 'xr', markersize = 2.0, label = f'ROM-FOM e: {(np.linalg.norm(fom-rom)/np.linalg.norm(fom)):.2E}')
            
            #### RBF ####
            rbf_data_filename = f"RBF_interpolated_data/{mu[0]}, {mu[1]}.npy"
            if os.path.exists(rbf_data_filename):
                fom = np.load(f'DataBase/Snapshots/{mu[0]}, {mu[1]}.npy')
                rbf = np.load(rbf_data_filename).T

                LaunchFakeSimulation(rbf, mu)
                rbf_skin_data_filename = f"RBF_Skin_Data/{mu[0]}, {mu[1]}.dat"
                x_rbf  = np.loadtxt(rbf_skin_data_filename, usecols=(0,))
                cp_rbf = np.loadtxt(rbf_skin_data_filename, usecols=(3,))
                fig = plt.plot(x_rbf, cp_rbf, '+g', markersize = 2.0, label = f'RBF-FOM e: {(np.linalg.norm(fom-rbf)/np.linalg.norm(rbf)):.2E}')

            if np.min(cp_fom) < cp_min:
                cp_min = np.min(cp_fom)
            if np.max(cp_fom) > cp_max:
                cp_max = np.max(cp_fom)
            fig = plt.title('Cp vs x')
            fig = plt.axis([-0.05,1.35,cp_max+0.1,cp_min-0.1])
            fig = plt.ylabel('Cp')
            fig = plt.xlabel('x')
            fig = plt.grid()
            fig = plt.legend()
            fig = plt.tight_layout()
            fig = plt.savefig(capture_name)
            fig = plt.close('all')

    if os.path.exists('mu_test.npy') and test:
        mu_test = np.load('mu_test.npy')

        for mu in mu_test:

            rbf_data_filename = f"RBF_interpolated_data/{mu[0]}, {mu[1]}.npy"
            if os.path.exists(rbf_data_filename):

                if not os.path.exists('Test_fom_rom_rbf_errors'):
                    os.mkdir('Test_fom_rom_rbf_errors')

                if not os.path.exists('RBF_Skin_Data'):
                    os.mkdir('RBF_Skin_Data')

                capture_name = f"Test_fom_rom_rbf_errors/Test_{mu[0]}, {mu[1]}.png"

                #### ONLINE CP PLOT
                ######################################################################
                cp_min = 0
                cp_max = 0
                fig = plt.figure()
                fig.set_figwidth(12.0)
                fig.set_figheight(8.0)

                #### FOM ####
                fom_skin_data_filename = f"DataBase/FOM_Skin_Data/{mu[0]}, {mu[1]}.dat"
                if os.path.exists(fom_skin_data_filename):
                    x_fom  = np.loadtxt(fom_skin_data_filename, usecols=(0,))
                    cp_fom = np.loadtxt(fom_skin_data_filename, usecols=(3,))
                    fig = plt.plot(x_fom, cp_fom, 'ob', markersize = 2.0, label = 'FOM')

                #### ROM ####
                rom_skin_data_filename = f"Skin_Data/ROM_Test{mu[0]}, {mu[1]}.dat"
                if os.path.exists(rom_skin_data_filename):
                    fom = np.load(f'DataBase/Snapshots/{mu[0]}, {mu[1]}.npy')
                    rom = np.load(f'RomDataBase/{mu[0]}, {mu[1]}.npy')
                    x_rom  = np.loadtxt(rom_skin_data_filename, usecols=(0,))
                    cp_rom = np.loadtxt(rom_skin_data_filename, usecols=(3,))
                    fig = plt.plot(x_rom, cp_rom, 'xr', markersize = 2.0, label = f'ROM-FOM e: {(np.linalg.norm(fom-rom)/np.linalg.norm(fom)):.2E}')
                
                #### RBF ####
                fom = np.load(f'DataBase/Snapshots/{mu[0]}, {mu[1]}.npy')
                rbf = np.load(rbf_data_filename).T

                LaunchFakeSimulation(rbf, mu)
                rbf_skin_data_filename = f"RBF_Skin_Data/{mu[0]}, {mu[1]}.dat"
                x_rbf  = np.loadtxt(rbf_skin_data_filename, usecols=(0,))
                cp_rbf = np.loadtxt(rbf_skin_data_filename, usecols=(3,))
                fig = plt.plot(x_rbf, cp_rbf, '+g', markersize = 2.0, label = f'RBF-FOM e: {(np.linalg.norm(fom-rbf)/np.linalg.norm(rbf)):.2E}')

                if np.min(cp_fom) < cp_min:
                    cp_min = np.min(cp_fom)
                if np.max(cp_fom) > cp_max:
                    cp_max = np.max(cp_fom)
                fig = plt.title('Cp vs x')
                fig = plt.axis([-0.05,1.35,cp_max+0.1,cp_min-0.1])
                fig = plt.ylabel('Cp')
                fig = plt.xlabel('x')
                fig = plt.grid()
                fig = plt.legend()
                fig = plt.tight_layout()
                fig = plt.savefig(capture_name)
                fig = plt.close('all')