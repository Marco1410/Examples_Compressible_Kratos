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
        parameters_copy = FakeProjectParameters(parameters.Clone(), mu)
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
            nametype = parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].GetString()
            if 'HHROM' in nametype:
                simulation_name = nametype.removeprefix(f"Results/HHROM_")
                skin_data_filename = f"HHROM_Skin_Data/{simulation_name}.dat"
            else:
                simulation_name = nametype.removeprefix(f"Results/RBF_")
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
# load parameters
#
def load_mu_parameters():
    if os.path.exists("mu_train.npy") and os.path.exists("mu_test.npy"):
        mu_train = np.load('mu_train.npy')
        mu_test = np.load('mu_test.npy')
        mu_train =  [mu.tolist() for mu in mu_train]
        mu_test =  [mu.tolist() for mu in mu_test]
    elif os.path.exists("mu_train.npy"):
        mu_train = np.load('mu_train.npy')
        mu_train =  [mu.tolist() for mu in mu_train]
        mu_test = []
    elif os.path.exists("mu_test.npy"):
        mu_test = np.load('mu_test.npy')
        mu_test =  [mu.tolist() for mu in mu_test]
        mu_train = []
    return mu_train, mu_test

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# save parameters
#
def FakeProjectParameters(parameters, mu=None):
    angle_of_attack = mu[0]
    mach_infinity   = mu[1]
    parameters["solver_settings"]["echo_level"].SetInt(0)
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(angle_of_attack)
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach_infinity)
    parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString(f'Results/{mu[2]}_{angle_of_attack}, {mach_infinity}')
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
# save parameters
#
def Plot_Cps(mu_list, capture_directory):
    #Disable logs
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

    for mu in mu_list:
        case_name = f'{mu[0]}, {mu[1]}'
        capture_filename   = f"{capture_directory}/{case_name}.png"

        #### CP PLOT
        ######################################################################
        cp_min = cp_max = cp_hrom = cp_hhrom = cp_fom = cp_rom = cp_rbf = 0
        fig = plt.figure()
        fig.set_figwidth(12.0)
        fig.set_figheight(8.0)

        if capture_directory == 'Validation':
            #### VALIDATION ######
            validation_skin_data_filename = f"../reference_data/flo36/{case_name}.dat"
            if os.path.exists(validation_skin_data_filename):
                x  = np.loadtxt(validation_skin_data_filename, usecols=(0,))
                cp_validation = np.loadtxt(validation_skin_data_filename, usecols=(1,))
                fig = plt.plot(x, cp_validation, 'r.-', markersize = 5.0, label = 'flo36')
            #### XFLR5 ######
            validation_skin_data_filename = f"../reference_data/flo36/xflr5_polars/{case_name}.dat"
            if os.path.exists(validation_skin_data_filename):
                x_xfoil  = np.loadtxt(validation_skin_data_filename, usecols=(0,))
                cp_xfoil = np.loadtxt(validation_skin_data_filename, usecols=(1,))
                fig = plt.plot(x_xfoil, cp_xfoil, '*', markersize = 5.0, label = 'XFLR5')
        #### FOM ######
        fom_skin_data_filename = f"FOM_Skin_Data/{case_name}.dat"
        if os.path.exists(fom_skin_data_filename):
            x_fom  = np.loadtxt(fom_skin_data_filename, usecols=(0,))
            y_fom  = np.loadtxt(fom_skin_data_filename, usecols=(1,))
            cp_fom = np.loadtxt(fom_skin_data_filename, usecols=(3,))
            if capture_directory == 'Validation':
                y_gtz = y_fom > 0
                y_lez = y_fom <= 0
                x_gtz  = x_fom[y_gtz] 
                cp_gtz = cp_fom[y_gtz]
                indices_ordenados = sorted(range(len(x_gtz)), key=lambda i: x_gtz[i])
                x_gtz = [x_gtz[i] for i in indices_ordenados]
                cp_gtz = [cp_gtz[i] for i in indices_ordenados]
                x_lez  = x_fom[y_lez] 
                cp_lez = cp_fom[y_lez]
                indices_ordenados = sorted(range(len(x_lez)), key=lambda i: x_lez[i])
                x_lez = [x_lez[i] for i in indices_ordenados]
                cp_lez = [cp_lez[i] for i in indices_ordenados]
                cp_val_gtz = np.interp(x_gtz, x[79:][::-1], cp_validation[79:][::-1])
                cp_val_lez = np.interp(x_lez, x[:78], cp_validation[:78])
                x_sim = np.concatenate((x_lez,x_gtz[::-1]))
                cp_sim = np.concatenate((cp_lez, cp_gtz[::-1]))
                cp_val = np.concatenate((cp_val_lez,cp_val_gtz[::-1]))
                error = (np.linalg.norm(cp_sim-cp_val)/np.linalg.norm(cp_val))
                fig = plt.plot(x_fom, cp_fom, 's', markersize = 5.0, label = f'FOM-Validation e: {error:.2E}')
            else:
                fig = plt.plot(x_fom, cp_fom, 's', markersize = 5.0, label = f'FOM')
        #### ROM ######
        rom_skin_data_filename = f"ROM_Skin_Data/{case_name}.dat"
        if os.path.exists(rom_skin_data_filename):
            x_rom  = np.loadtxt(rom_skin_data_filename, usecols=(0,))
            cp_rom = np.loadtxt(rom_skin_data_filename, usecols=(3,))
            fig = plt.plot(x_rom, cp_rom, 'o', markersize = 3.0, label = f'ROM-FOM e: {(np.linalg.norm(cp_fom-cp_rom)/np.linalg.norm(cp_fom)):.2E}')
        #### HROM ######
        hrom_skin_data_filename = f"HROM_Skin_Data/{case_name}.dat"
        if os.path.exists(hrom_skin_data_filename):
            x_hrom  = np.loadtxt(hrom_skin_data_filename, usecols=(0,))
            cp_hrom = np.loadtxt(hrom_skin_data_filename, usecols=(3,))
            fig = plt.plot(x_hrom, cp_hrom, '+', markersize = 6.0, label = f'HROM-FOM e: {(np.linalg.norm(cp_fom-cp_hrom)/np.linalg.norm(cp_fom)):.2E}')
        #### HHROM ######
        hhrom_name = f"HHROM_Snapshots/{case_name}.npy"
        if os.path.exists(hhrom_name):
            hhrom = np.load(hhrom_name)
            LaunchFakeSimulation(hhrom, [mu[0], mu[1], 'HHROM'])
            hhrom_skin_data_filename = f"HHROM_Skin_Data/{case_name}.dat"
            x_hhrom  = np.loadtxt(hhrom_skin_data_filename, usecols=(0,))
            cp_hhrom = np.loadtxt(hhrom_skin_data_filename, usecols=(3,))
            fig = plt.plot(x_hhrom, cp_hhrom, 'x', markersize = 6.0, label = f'HHROM-FOM e: {(np.linalg.norm(cp_fom-cp_hhrom)/np.linalg.norm(cp_fom)):.2E}')
        #### RBF ####
        rbf_name = f"RBF_Snapshots/{case_name}.npy"
        if os.path.exists(rbf_name):
            rbf = np.load(rbf_name).T
            LaunchFakeSimulation(rbf, [mu[0], mu[1], 'RBF'])
            rbf_skin_data_filename = f"RBF_Skin_Data/{case_name}.dat"
            x_rbf  = np.loadtxt(rbf_skin_data_filename, usecols=(0,))
            cp_rbf = np.loadtxt(rbf_skin_data_filename, usecols=(3,))
            fig = plt.plot(x_rbf, cp_rbf, '.', markersize = 2.0, label = f'RBF-FOM e: {(np.linalg.norm(cp_fom-cp_rbf)/np.linalg.norm(cp_fom)):.2E}')
        
        cp_min = np.min([np.min(cp_hhrom), np.min(cp_hrom), np.min(cp_fom), np.min(cp_rom), np.min(cp_rbf)])
        cp_max = np.max([np.max(cp_hhrom), np.max(cp_hrom), np.max(cp_fom), np.max(cp_rom), np.max(cp_rbf)])
        
        fig = plt.title('Cp vs x')
        fig = plt.axis([-0.05,1.05,cp_max+0.1,cp_min-0.1])
        fig = plt.ylabel('Cp')
        fig = plt.xlabel('x')
        fig = plt.grid()
        fig = plt.legend()
        fig = plt.tight_layout()
        fig = plt.savefig(capture_filename, dpi=400)
        fig = plt.close('all')


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


if __name__ == "__main__":

    mu_train, mu_test = load_mu_parameters()

    Plot_Cps(mu_train, 'Train_Captures')
    Plot_Cps(mu_test, 'Test_Captures')

    mu = []
    mu.append([1.0,0.72])
    mu.append([1.0,0.73])
    mu.append([1.0,0.75])
    mu.append([2.0,0.75])

    Plot_Cps(mu, 'Validation')