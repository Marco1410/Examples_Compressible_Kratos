import os
import importlib
import numpy as np
import KratosMultiphysics
import matplotlib.pyplot as plt
import KratosMultiphysics.kratos_utilities
from KratosMultiphysics.gid_output_process import GiDOutputProcess
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition


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

            simulation_name = nametype.removeprefix(f"Reconstruction_Test_Captures/Results/RB_")
            skin_data_filename = f"Reconstruction_Test_Captures/RB_Skin_Data/{simulation_name}.dat"

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
def FakeProjectParameters(parameters, mu=None):
    angle_of_attack = mu[0]
    mach_infinity   = mu[1]
    parameters["solver_settings"]["echo_level"].SetInt(0)
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(angle_of_attack)
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach_infinity)
    parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString(f'Reconstruction_Test_Captures/Results/{mu[2]}_{angle_of_attack}, {mach_infinity}')
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
# Plot Cps 
#
def Plot_Cps(mu_list, capture_directory):
    #Disable logs
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

    for mu in mu_list:
        case_name = f'{mu[0]}, {mu[1]}'
        capture_filename   = f"{capture_directory}/Plots/{case_name}.png"

        #### CP PLOT
        ######################################################################
        cp_min = cp_max = cp_fom = cp_rebuild = 0
        fig = plt.figure()
        fig.set_figwidth(12.0)
        fig.set_figheight(8.0)

        #### FOM ######
        fom_skin_data_filename = f"FOM_Skin_Data/{case_name}.dat"
        fom_data_filename      = f"FOM_Snapshots/{case_name}.npy"
        if os.path.exists(fom_skin_data_filename) and os.path.exists(fom_data_filename):
            fom    = np.load(fom_data_filename)
            x_fom  = np.loadtxt(fom_skin_data_filename, usecols=(0,))
            cp_fom = np.loadtxt(fom_skin_data_filename, usecols=(3,))
            fig = plt.plot(x_fom, cp_fom, 's', markersize = 5.0, label = f'FOM')

        #### REBUILD ####
        rb_skin_data_filename = f"{capture_directory}/RB_Skin_Data/{case_name}.dat"
        rb_data_filename      = f"{capture_directory}/RB_Snapshots/{case_name}.npy"
        if os.path.exists(rb_data_filename):
            rebuild = np.load(rb_data_filename)
            LaunchFakeSimulation(rebuild, [mu[0], mu[1], 'RB'])
            x_rebuild  = np.loadtxt(rb_skin_data_filename, usecols=(0,))
            cp_rebuild = np.loadtxt(rb_skin_data_filename, usecols=(3,))

            error_t  = np.linalg.norm(fom-rebuild)/np.linalg.norm(fom)
            error_t  = np.linalg.norm(fom-rebuild)/np.linalg.norm(fom)

            error_cp = np.linalg.norm(cp_fom-cp_rebuild)/np.linalg.norm(cp_fom)
            error_cp = np.linalg.norm(cp_fom-cp_rebuild)/np.linalg.norm(cp_fom)

            fig = plt.plot(x_rebuild, cp_rebuild, '.', markersize = 2.0, label = f'RB-FOM e_t : {error_t:.2E} e_cp: {error_cp:.2E}')
        
        cp_min = np.min([np.min(cp_fom), np.min(cp_rebuild)])
        cp_max = np.max([np.max(cp_fom), np.max(cp_rebuild)])

        fig = plt.title('Cp vs x')
        fig = plt.axis([-0.05,1.05,cp_max+0.1,cp_min-0.1])
        fig = plt.axis()
        fig = plt.ylabel('Cp')
        fig = plt.xlabel('x')
        fig = plt.grid()
        fig = plt.legend()
        fig = plt.tight_layout()
        # fig = plt.show()
        fig = plt.savefig(capture_filename, dpi=400)
        fig = plt.close('all')


def DEIM(Basis):
    #find first point
    U = Basis[:,0].reshape(Basis.shape[0], -1)
    z = np.zeros(U.shape)
    P = z.copy()
    indexes = np.argmax( np.abs(Basis[:,0]) )
    P[indexes] = 1

    #find next points
    for i in range(1,Basis.shape[1]):
        c = np.linalg.solve(P.T @ U , P.T @ Basis[:,i] )
        residual = Basis[:,i] - U @ c
        U = np.c_[U,Basis[:,i]]
        index_i = np.argmax( np.abs(residual) )
        indexes = np.r_[indexes, index_i]
        P = np.c_[P, z]; P[index_i,i] = 1

    return indexes

if '__main__':

    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Reconstruction_Test_Captures')

    folder_names = ["Reconstruction_Test_Captures", "Reconstruction_Test_Captures/Plots",
                    "Reconstruction_Test_Captures/RB_Skin_Data", "Reconstruction_Test_Captures/RB_Snapshots"]
    
    for name in folder_names:
        if not os.path.exists(name):
            os.mkdir(name)

    mu_train, mu_test = load_mu_parameters()

    snapshots = []
    for mu in mu_train:
        file = f'{mu[0]}, {mu[1]}'
        snapshots.append(np.load(f'FOM_Snapshots/{file}.npy'))
        # snapshots.append(np.load(f'Q_errors/{file}.npy'))
    snapshots = np.block(snapshots)

    # Calculate the randomized SVD of the snapshots matrix
    svd_truncation_tolerance = 1e-6
    phi,_,_,_= RandomizedSingularValueDecomposition().Calculate(snapshots, svd_truncation_tolerance)

    # phi = np.load('rom_data/RightBasisMatrix.npy')

    rebuild = phi @ (phi.T @ snapshots)

    for i, mu in enumerate(mu_train):
        filename = f"Reconstruction_Test_Captures/RB_Snapshots/{mu[0]}, {mu[1]}"
        np.save(filename, rebuild[:,i].reshape(-1,1))

    Plot_Cps(mu_train, 'Reconstruction_Test_Captures')

    #using DEIM to select the best points
    DEIM_points = DEIM(phi)
    u_DEIM = phi[DEIM_points, :]

    np.save('rom_data/RightBasisMatrix', u_DEIM)

    S_DEIM = snapshots[DEIM_points, :]
    S_reconstructed_DEIM = phi @ (np.linalg.pinv(u_DEIM) @ S_DEIM)

    print('Number of modes taken:', phi.shape[1])
    print('Number of snapshots:', snapshots.shape[1])
    print(f'Total error: {np.linalg.norm(snapshots-rebuild)/np.linalg.norm(snapshots):.2E}')
    print('DEIM_points:', DEIM_points)
    print(f'Total DEIM error: {np.linalg.norm(snapshots-S_reconstructed_DEIM)/np.linalg.norm(snapshots):.2E}')