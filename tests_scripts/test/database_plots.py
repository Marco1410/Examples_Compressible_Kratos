import matplotlib.pyplot as plt
from parameters_manager import *  
from database_rbf import *
import importlib
import numpy as np
import KratosMultiphysics
import KratosMultiphysics.kratos_utilities
from KratosMultiphysics.RomApplication.rom_database import RomDatabase
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
from KratosMultiphysics.gid_output_process import GiDOutputProcess
from KratosMultiphysics.vtk_output_process import VtkOutputProcess




def plot_parameters_with_cp(target_directory='output_plots', mu_id=0, svd_tol=1e-12, weight=1.0, alpha=1.0, beta=1.0,indices=[0,0,0,0]):
    os.makedirs(target_directory, exist_ok=True)
    mu_test  = [[mu[0],mu[1]] for mu in load_mu_parameters('mu_test')]   
    mu_validation = [[mu[0],mu[1]] for mu in load_mu_parameters('mu_validation')]  
    mu_train = [[mu[0],mu[1]] for mu in load_mu_parameters(f'mu_train_{mu_id}_{alpha}_{beta}')]  
    # Initialize database
    data_base = initialize_database(mu_id=mu_id,weight=weight,svd_tol=svd_tol,alpha=alpha,beta=beta) 
    mach_train = [mu[1] for mu in mu_train]
    alpha_train = [mu[0] for mu in mu_train]
    mach_test = [mu[1] for mu in mu_test]
    alpha_test = [mu[0] for mu in mu_test]
    mach_validation = [mu[1] for mu in mu_validation]
    alpha_validation = [mu[0] for mu in mu_validation]
    fig = plt.figure(figsize=[9,8])
    central_plot_position = [0.32, 0.32, 0.4, 0.4]  # [left, bottom, width, height]
    ax = fig.add_axes(central_plot_position)  
    ax.plot(mach_train     , alpha_train     , 'bs', label='Train'     , markersize=4)
    ax.plot(mach_test      , alpha_test      , 'rx', label='Test'      , markersize=7)
    ax.plot(mach_validation, alpha_validation, 'g*', label='Validation', markersize=10)
    ax.set_xlabel('Mach')
    ax.set_ylabel('Alpha')
    ax.grid(True, linestyle='--', alpha=0.8)
    positions = [
        [0.08, 0.75, 0.20, 0.20],
        [0.76, 0.75, 0.20, 0.20],
        [0.08, 0.06, 0.20, 0.20],
        [0.76, 0.06, 0.20, 0.20]]
    for idx, pos in zip(indices,positions):
        mu = [alpha_train[idx],mach_train[idx],weight,beta,mu_id,svd_tol]
        in_database, hash_cp_data = data_base.check_if_in_database("QoI_FOM", mu)
        if in_database:
            cp_data = data_base.get_single_numpy_from_database(f'{hash_cp_data}Cp_data')
            sub_ax = fig.add_axes(pos)
            x  = cp_data[:,0]
            cp = cp_data[:,3]
            sub_ax.plot( x, -cp, '.', markersize = 1.0)
            sub_ax.set_title(f'Angle={np.round(alpha_train[idx],3)}, Mach={np.round(mach_train[idx],3)}', size=13, pad=8)
            sub_ax.set_ylim(-1.5,1.5)
            sub_ax.set_xlabel('x', size=13)
            sub_ax.set_ylabel('$-C_p$', size=13)
            arrowprops = dict(arrowstyle="->", color='grey', lw=1.5)
            ax.annotate('',
                        xy=(mach_train[idx], alpha_train[idx]), 
                        xycoords='data',
                        xytext=(pos[0] + pos[2]/2, pos[1] + pos[3]/2), 
                        textcoords='figure fraction',
                        arrowprops=arrowprops)
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, fontsize=12, bbox_to_anchor=(0.515, 0.79))
    plt.savefig(f'{target_directory}/Mu_parameters_id_{mu_id}.png', dpi=400)
    plt.close('all')


def plot_parameters(target_directory='output_plots', mu_ids=[0], alpha_values=[1.0]):
    os.makedirs(target_directory, exist_ok=True)
    mu_test  = [[mu[0], mu[1]] for mu in load_mu_parameters('mu_test')]   
    mu_validation = [[mu[0], mu[1]] for mu in load_mu_parameters('mu_validation')] 
    mach_test = [mu[1] for mu in mu_test]
    alpha_test = [mu[0] for mu in mu_test]
    mach_validation = [mu[1] for mu in mu_validation]
    alpha_validation = [mu[0] for mu in mu_validation]
    # fig, axes = plt.subplots(nrows=1, ncols=len(mu_ids), figsize=(10, 5), sharex=True, sharey=True)
    fig, axes = plt.subplots(nrows=1, ncols=len(mu_ids), sharex=True, sharey=True)
    handles, labels = None, None 
    for id, (mu_id, alpha) in enumerate(zip(mu_ids, alpha_values)):
        mu_train = [[mu[0], mu[1]] for mu in load_mu_parameters(f'mu_train_{mu_id}_{alpha}_{alpha}')]  
        mach_train = [mu[1] for mu in mu_train]
        alpha_train = [mu[0] for mu in mu_train]
        if len(mu_ids) == 1:
            p1 = axes.plot(mach_train, alpha_train, 'bs', label='Train', markersize=4)
            p2 = axes.plot(mach_test, alpha_test, 'rx', label='Test', markersize=7)
            p3 = axes.plot(mach_validation, alpha_validation, 'g*', label='Validation', markersize=10)
            axes.grid(True, linestyle='--', alpha=0.8)
            axes.set_title(f'Density constant = {alpha}', size=13, pad=8)
            if handles is None and labels is None:
                handles, labels = axes.get_legend_handles_labels()
        else:
            p1 = axes[id].plot(mach_train, alpha_train, 'bs', label='Train', markersize=4)
            p2 = axes[id].plot(mach_test, alpha_test, 'rx', label='Test', markersize=7)
            p3 = axes[id].plot(mach_validation, alpha_validation, 'g*', label='Validation', markersize=10)
            axes[id].grid(True, linestyle='--', alpha=0.8)
            axes[id].set_title(f'Density constant = {alpha}', size=13, pad=8)
            if handles is None and labels is None:
                handles, labels = axes[id].get_legend_handles_labels()
    fig.supxlabel('Mach', fontsize=14, y=0.04)
    fig.supylabel('Alpha', fontsize=14, x=0.015)
    fig.suptitle(f'Parameter space m = {len(mu_train)} - mu {mu_id}', size=16)
    fig.legend(handles, labels, loc="upper center", ncol=3, fontsize=12, bbox_to_anchor=(0.5, 0.92))
    plt.tight_layout(rect=[0, 0, 1, 0.92])
    plt.savefig(f'{target_directory}/Mu_parameters.png', dpi=400)
    plt.close('all')


def plot_decay(target_directory='output_plots', weight=1.0, mu_ids=[0], svd_tols=[1e-12], alpha_values=[1.0]):
    os.makedirs(target_directory, exist_ok=True) 
    for (mu_id, svd_tol, alpha) in zip(mu_ids, svd_tols, alpha_values):
        # Initialize database
        data_base = initialize_database(mu_id=mu_id,weight=weight,svd_tol=svd_tol,alpha=alpha,beta=alpha) 
        snapshots = data_base.get_snapshots_matrix_from_database(load_mu_parameters(f'mu_train_{mu_id}_{alpha}_{alpha}'), table_name='FOM')
        _,singular_values,_,_ = RandomizedSingularValueDecomposition(COMPUTE_V=False).Calculate(snapshots,svd_tol)
        singular_values = np.sort(singular_values)[::-1]
        acumulated_sum = np.cumsum(singular_values)
        normalized_acumulated_sum = acumulated_sum / acumulated_sum[-1]
        decay = 1 - normalized_acumulated_sum
        # plt.figure(figsize=(10, 6))
        plt.figure()
        plt.semilogy(range(len(singular_values)), decay, markersize=5, linestyle='', color='b', label=f'Mu_id: {mu_id}')
    plt.title('Singular values decay')
    plt.xlabel('Singular value index')
    plt.ylabel('Decay')
    plt.legend()
    plt.grid(True, which="both", ls="--")
    plt.savefig(f'{target_directory}/Singular_values_decay_mu_train_{mu_id}_{alpha}_{alpha}.png', bbox_inches='tight')
    # plt.show()
    plt.close()


def plot_fom_rom_x_validation(target_directory='output_plots',
                            mu_id   = 0    ,  
                            svd_tol = 1e-12,  
                            weight  = 1.0  ,  
                            alpha   = 1.0  , 
                            beta    = 1.0  ):
    os.makedirs(target_directory, exist_ok=True)

    mu_validation = [[mu[0],mu[1]] for mu in load_mu_parameters('mu_validation')]  

    # Initialize database
    with open("rom_manager_parameters.json", 'r') as parameter_file:
        general_rom_manager_parameters = KratosMultiphysics.Parameters(parameter_file.read())
    general_rom_manager_parameters["ROM"]["svd_truncation_tolerance"].SetDouble(svd_tol)
    data_base = RomDatabase(general_rom_manager_parameters, mu_names=['Alpha','Mach','Weight','Beta','MuId','SVDTolerance'])

    alpha_values = np.array([mu[0] for mu in mu_validation])
    mach_values  = np.array([mu[1] for mu in mu_validation])

    # Necessary to search in the database
    mu = [alpha_values[0],mach_values[0],weight,beta,mu_id,svd_tol]
    in_database, hash_cp_data = data_base.check_if_in_database("QoI_FOM", mu)

    if in_database:
        cp_data = data_base.get_single_numpy_from_database(f'{hash_cp_data}Cp_data')
        cp = cp_data[:,3]

    Cp_FOM = np.loadtxt('FOM_0.45, 0.73.dat', usecols=(3,))
    Cp_ROM = np.loadtxt('ROM_0.45, 0.73.dat', usecols=(3,))

    plt.figure(figsize=(6, 5))
    plt.scatter(Cp_FOM, Cp_ROM, c='w', alpha=0.5, edgecolors='r', label="Data")
    val_max = np.max(Cp_FOM)
    val_min = np.min(Cp_FOM)
    plt.plot([val_min, val_max], [val_min, val_max], 'k--', linewidth=1, label="y=x")
    plt.xlabel("$C_p$ FOM")
    plt.ylabel("$C_p$ ROM")
    plt.title("ROM vs FOM")
    plt.grid()
    plt.legend()
    # plt.show()
    plt.savefig(f'{target_directory}/Cp_fom_rom_x_validation.png', dpi=400)
    plt.close('all')


##################################################################################################################################
def LaunchFakeSimulation(data_set, mode, path):
    with open('ProjectParameters.json','r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    model = KratosMultiphysics.Model()
    parameters_copy = FakeProjectParameters(parameters.Clone(), mode, path)
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
def FakeProjectParameters(parameters, mode, path):
    parameters["solver_settings"]["echo_level"].SetInt(0)
    parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString(f'{path}/mode_{mode}')
    parameters["output_processes"]["vtk_output"][0]["Parameters"]["output_path"].SetString(f'{path}/mode_{mode}')
    return parameters
def _GetAnalysisStageClass(parameters):
    analysis_stage_module_name = parameters["analysis_stage"].GetString()
    analysis_stage_class_name = analysis_stage_module_name.split('.')[-1]
    analysis_stage_class_name = ''.join(x.title() for x in analysis_stage_class_name.split('_'))
    analysis_stage_module = importlib.import_module(analysis_stage_module_name)
    analysis_stage_class = getattr(analysis_stage_module, analysis_stage_class_name)
    return analysis_stage_class
##################################################################################################################################


def plot_modes_visualization(modes_to_plot=[0,1,2,3,4], target_directory='output_plots', 
                             mu_id=0, svd_tol=1e-12, weight=1.0, alpha=1.0, beta=1.0):
    os.makedirs(target_directory, exist_ok=True) 
    # Initialize database
    data_base = initialize_database(mu_id=mu_id,weight=weight,svd_tol=svd_tol,alpha=alpha,beta=beta) 
    mu_train = [[mu[0],mu[1]] for mu in load_mu_parameters(f'mu_train_{mu_id}_{alpha}_{beta}')]  
    in_database, hash_basis = data_base.check_if_in_database("RightBasis", mu_train)
    if not in_database:
        print(f'No RightBasis data for mu_train_{mu_id}_{alpha}_{beta}.')
    else:
        phi = data_base.get_single_numpy_from_database(hash_basis)
        for mode in modes_to_plot:
            mode_values = phi[:,mode]
            LaunchFakeSimulation(mode_values, mode, target_directory)



# def plot_Cps(mu, capture_directory, fom_skin_cp, rbf_skin_cp_from_variables, rbf_skin_cp_from_full_cp, rbf_skin_cp):
#     os.makedirs(capture_directory, exist_ok=True)
#     case_name = f"{mu[0]}, {mu[1]}"
#     capture_filename = f"{capture_directory}/{case_name}.png"
#     plt.figure(figsize=(9, 6))
#     cp_validation = 0
#     # VALIDATION
#     if 'validation' in capture_directory:
#         validation_filename = f"reference_data/flo36/{case_name}.dat"
#         if os.path.exists(validation_filename):
#             x = np.loadtxt(validation_filename, usecols=(0,))
#             cp_validation = np.loadtxt(validation_filename, usecols=(1,))
#             plt.plot(x, cp_validation, '*-', markersize=10.0, label='flo36')
#     # FOM
#     x_fom, y_fom, cp_fom = fom_skin_cp[:, 0].reshape(-1,1), fom_skin_cp[:, 1].reshape(-1,1), fom_skin_cp[:, 3].reshape(-1,1)
#     if 'validation' in capture_directory:
#         mask_gtz = y_fom > 0
#         mask_lez = ~mask_gtz
#         x_sorted_gtz = x_fom[mask_gtz][np.argsort(x_fom[mask_gtz])]
#         cp_sorted_gtz = cp_fom[mask_gtz][np.argsort(x_fom[mask_gtz])]
#         x_sorted_lez = x_fom[mask_lez][np.argsort(x_fom[mask_lez])]
#         cp_sorted_lez = cp_fom[mask_lez][np.argsort(x_fom[mask_lez])]
#         cp_val_gtz = np.interp(x_sorted_gtz, x[79:][::-1], cp_validation[79:][::-1])
#         cp_val_lez = np.interp(x_sorted_lez, x[:78], cp_validation[:78])
#         cp_sim = np.concatenate((cp_sorted_lez, cp_sorted_gtz[::-1]))
#         cp_val = np.concatenate((cp_val_lez, cp_val_gtz[::-1]))
#         error = np.linalg.norm(cp_sim - cp_val) / np.linalg.norm(cp_val)
#         plt.plot(x_fom, cp_fom, 's', markersize=3.5, label=f'FOM-Validation e: {error:.2E}')
#     else:
#         plt.plot(x_fom, cp_fom, 's', markersize=3.5, label='FOM')
#     # RBF Comparisons
#     x_rbf, cp_rbf_from_variables = rbf_skin_cp_from_variables[:, 0].reshape(-1,1), rbf_skin_cp_from_variables[:, 3].reshape(-1,1)
#     plt.plot(x_rbf, cp_rbf_from_variables, '.', markersize=8.0, label=f'RBF_pot_vars-FOM e: {(np.linalg.norm(cp_fom - cp_rbf_from_variables) / np.linalg.norm(cp_fom)):.2E}')
#     plt.plot(x_rbf, rbf_skin_cp_from_full_cp, '.', markersize=6.0, label=f'RBF_full_cp-FOM e: {(np.linalg.norm(cp_fom - rbf_skin_cp_from_full_cp) / np.linalg.norm(cp_fom)):.2E}')
#     plt.plot(x_rbf, rbf_skin_cp, '.', markersize=2.0, label=f'RBF_skin_cp-FOM e: {(np.linalg.norm(cp_fom - rbf_skin_cp) / np.linalg.norm(cp_fom)):.2E}')
#     cp_min = np.nanmin([np.min(cp_fom), np.min(cp_rbf_from_variables), np.min(rbf_skin_cp_from_full_cp), np.min(rbf_skin_cp), np.min(cp_validation)])
#     cp_max = np.nanmax([np.max(cp_fom), np.max(cp_rbf_from_variables), np.max(rbf_skin_cp_from_full_cp), np.max(rbf_skin_cp), np.max(cp_validation)])
#     plt.axis([-0.05,1.05,cp_max+0.1,cp_min-0.1])
#     plt.title('Cp vs x')
#     plt.xlabel('x')
#     plt.ylabel('-Cp')
#     plt.legend()
#     plt.grid()
#     plt.tight_layout()
#     plt.savefig(capture_filename, dpi=200)
#     plt.close('all')


# def plot_Cps_2d(mu, capture_directory, case_types=['FOM'], rbf_path='output_plots/RBF/potential_variables'):
#     os.makedirs(capture_directory, exist_ok=True)
#     case_name = f"{mu[0]}, {mu[1]}.png"
#     simulation_type = next((s for s in ['Fit', 'Test', 'Run'] if s in capture_directory), 'Unknown')
#     for case in case_types:
#         path = f"{rbf_path}/Results/vtk_output_{mu[0]}, {mu[1]}" if 'RBF' in case else f"Results/vtk_output_{case}_{simulation_type}{', '.join(map(str, mu))}"
#         if os.path.exists(path):
#             vtk_filename = get_latest_vtk_file(path)
#             if os.path.exists(vtk_filename):
#                 mesh = pv.read(vtk_filename).extract_geometry()
#                 points, cells = mesh.points, mesh.faces.reshape(-1, 4)[:, 1:4]
#                 if 'full_cp' in rbf_path:
#                     pressure_coeff = mesh.point_data['REACTION_WATER_PRESSURE']
#                 else:
#                     pressure_coeff = mesh.point_data['PRESSURE_COEFFICIENT']
#                 fig, ax = plt.subplots(figsize=(10, 8))
#                 triang = tri.Triangulation(points[:, 0], points[:, 1], cells)
#                 cmap = cm.get_cmap('viridis') # viridis coolwarm
#                 contour = ax.tricontour(triang, pressure_coeff, levels=15, cmap='coolwarm', linewidths=0.8)
#                 ax.tricontourf(triang, pressure_coeff, levels=100, cmap=cmap, alpha=1.0)
#                 cbar = plt.colorbar(cm.ScalarMappable(norm=plt.Normalize(vmin=pressure_coeff.min(), vmax=pressure_coeff.max()), cmap=cmap), ax=ax, label='PRESSURE COEFFICIENT')
#                 cbar.ax.tick_params(labelsize=10)
#                 ax.set_xlabel('X')
#                 ax.set_ylabel('Y')
#                 ax.set_title(f'{case} Pressure Coefficient Distribution - Alpha={mu[0]:.3f}, Mach={mu[1]:.3f}', size=13, pad=10)
#                 ax.set_xlim(-0.5, 1.5)
#                 ax.set_ylim(-0.75, 1.0)
#                 plt.savefig(f'{capture_directory}/2D_{case}_{case_name}', dpi=200)
#                 plt.close('all')


def plot_captures(mu_list, target_directory, case, **params):
    data_base = initialize_database(mu_id=mu_id,weight=weight,svd_tol=svd_tol,alpha=alpha,beta=beta)  
    for mu in mu_list:
        mu_extended = [mu[0], mu[1], params['weight'], params['beta'], params['mu_id'], params['svd_tol'], params['solver_strategy']]
        if not (in_database := data_base.check_if_in_database("QoI_FOM", mu_extended))[0]:
            print(f'No Cp data for {mu}.')
            continue
        fom_skin_cp                = data_base.get_single_numpy_from_database(f'{in_database[1]}Cp_data')
        rbf_skin_cp_from_variables = np.load(f'{target_directory}/potential_variables/RBF_Skin_Data_skin_cp_from_potential_variables/{mu[0]}, {mu[1]}.npy')
        rbf_skin_cp_from_full_cp   = np.load(f'{target_directory}/full_cp/RBF_Skin_Data_skin_cp_from_full_cp/{mu[0]}, {mu[1]}.npy').reshape(-1,1)
        rbf_skin_cp                = np.load(f'{target_directory}/skin_cp/RBF_Skin_Data_skin_cp/{mu[0]}, {mu[1]}.npy').reshape(-1,1)
        plot_Cps(mu_extended, f'{target_directory}/{case}', fom_skin_cp, rbf_skin_cp_from_variables, rbf_skin_cp_from_full_cp, rbf_skin_cp)

        plot_Cps_2d(mu_extended, f'{target_directory}/{case}', case_types=['FOM'])
        plot_Cps_2d(mu_extended, f'{target_directory}/{case}', case_types=['RBF_potential_variables'], rbf_path=f'{target_directory}/potential_variables')
        plot_Cps_2d(mu_extended, f'{target_directory}/{case}', case_types=['RBF_full_cp'], rbf_path=f'{target_directory}/full_cp')





if __name__ == "__main__":
    #Disable logs
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

    target_directory = 'output_plots'
    os.makedirs(target_directory, exist_ok=True)

    plot_parameters_with_cps = False
    plot_parameters_sets     = False
    plot_modes               = False
    plot_decays              = False
    plot_cps_captures        = True

    # To plot parameters with Cps
    if plot_parameters_with_cps:
        mu_id   = 0      
        weight  = 1.0    
        svd_tol = 1e-12  
        alpha   = 1.0   
        beta    = 1.0   
        plot_parameters_with_cp(target_directory=target_directory, mu_id=mu_id, weight=weight, svd_tol=svd_tol, alpha=alpha, beta=beta, indices=[0,2,3,1])

    # Plot mu parameters
    if plot_parameters_sets:
        mu_ids          = [0]      
        alpha_values    = [1.0]   
        plot_parameters(target_directory=target_directory, mu_ids=mu_ids, alpha_values=alpha_values)

    # To plot modes
    if plot_modes:
        mu_id           = 0      
        weight          = 1.0    
        svd_tol         = 1e-12  
        alpha           = 1.0   
        beta            = 1.0 
        plot_modes_visualization(modes_to_plot=[0,1,2,3,4], target_directory=f'{target_directory}/Modes', 
                                mu_id=mu_id, weight=weight, svd_tol=svd_tol, alpha=alpha, beta=beta)

    # To plot svd decay
    if plot_decays:
        weight = 1.0
        mu_ids = [0]
        svd_tols = [1e-12]
        alpha_values = [1.0]
        plot_decay(target_directory=target_directory, weight=weight, mu_ids=mu_ids, svd_tols=svd_tols, alpha_values=alpha_values)

    if plot_cps_captures:
        directories = ['Training', 'Test', 'Validation']
        for directory in directories:
            KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting(f'{target_directory}/{directory}')
            os.makedirs(f'{target_directory}/{directory}', exist_ok=True)

        mu_train      = [[mu[0],mu[1]] for mu in load_mu_parameters(f'mu_train_{mu_id}_{alpha}_{beta}')]  
        mu_test       = [[mu[0],mu[1]] for mu in load_mu_parameters('mu_test')]   
        mu_validation = [[mu[0],mu[1]] for mu in load_mu_parameters('mu_validation')]  

        plot_captures(mu_list=mu_train, target_directory=target_directory,
                            case='Fit', mu_id=mu_id, weight=weight, svd_tol=svd_tol, alpha=alpha, beta=beta, solver_strategy=solver_strategy)

        # plot_Cps_distributions(mu_list=mu_test, target_directory=target_directory,
        #                     case='Test', mu_id=mu_id, weight=weight, svd_tol=svd_tol, alpha=alpha, beta=beta, solver_strategy=solver_strategy)

        # plot_Cps_distributions(mu_list=mu_validation, target_directory=target_directory,
        #                     case='Run_validation', mu_id=mu_id, weight=weight, svd_tol=svd_tol, alpha=alpha, beta=beta, solver_strategy=solver_strategy)
        
        # plot_fom_rom_x_validation(target_directory=f'{target_directory}/Dataset',
        #                                mu_id   = mu_id  , 
        #                                weight  = weight , 
        #                                svd_tol = svd_tol, 
        #                                alpha   = alpha  ,
        #                                beta    = beta   )