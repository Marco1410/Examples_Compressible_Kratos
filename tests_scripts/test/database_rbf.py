import os
import re
import importlib
import numpy as np
import pandas as pd
import pyvista as pv 
import seaborn as sns
import scipy.interpolate as interp
import matplotlib.cm as cm
import matplotlib.tri as tri
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
from ezyrb import ReducedOrderModel as ROM
from ezyrb import RBF, POD, Database
from parameters_manager import * 
import KratosMultiphysics
import KratosMultiphysics.kratos_utilities
from KratosMultiphysics.RomApplication.rom_database import RomDatabase
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
from KratosMultiphysics.gid_output_process import GiDOutputProcess
from KratosMultiphysics.vtk_output_process import VtkOutputProcess




def initialize_database(mu_id=0, weight=1, svd_tol=1e-12, alpha=1.0, beta=1.0):
    with open("rom_manager_parameters.json", 'r') as parameter_file:
        general_rom_manager_parameters = KratosMultiphysics.Parameters(parameter_file.read())
    general_rom_manager_parameters["ROM"]["svd_truncation_tolerance"].SetDouble(svd_tol)
    return RomDatabase(general_rom_manager_parameters, mu_names=['Alpha', 'Mach', 'Weight', 'Beta', 'MuId', 'SVDTolerance', 'RomSolverStrategy'])


def initialize_rbf(parameters, snapshots):
    #### RBF TRAINING
    return ROM(Database(parameters, snapshots.T), POD(), RBF()).fit()


def plot_L2_error_heat_map_interpolated(plotname='', mu_list=[], errors=None):
    alpha_values = np.array([mu[0] for mu in mu_list])
    mach_values  = np.array([mu[1] for mu in mu_list])
    errors = np.array(errors)
    alpha_grid = np.linspace(min(alpha_values), max(alpha_values), 100)
    mach_grid = np.linspace(min(mach_values), max(mach_values), 100)
    Mach, Alpha = np.meshgrid(mach_grid, alpha_grid)
    Error_grid = interp.griddata((mach_values, alpha_values), errors, (Mach, Alpha), method='cubic')
    # plt.figure(figsize=(8, 6))
    plt.figure()
    contour = plt.contourf(Mach, Alpha, Error_grid, levels=50, cmap="viridis") 
    plt.scatter(mach_values, alpha_values, edgecolors="k", s=50, alpha=0.75)
    plt.colorbar(contour, label=r"$L^2$ Error: $\frac{\| u - u_{\mathrm{ref}} \|_2}{\| u_{\mathrm{ref}} \|_2}$")
    plt.xlabel("Mach")
    plt.ylabel("Alpha")
    plt.title("$L^2$ error heat map")
    plt.grid(True, linestyle="--", alpha=0.5)
    # plt.show()
    plt.savefig(f'{plotname}.png', dpi=400)
    plt.close('all')


def plot_L2_error_heat_map_discreted(plotname='', mu_list=[], errors=None):
    alpha_values = np.array([mu[0] for mu in mu_list])
    mach_values  = np.array([mu[1] for mu in mu_list])
    errors = np.array(errors)
    # plt.figure(figsize=(8, 6))
    plt.figure()
    sc = plt.scatter(mach_values, alpha_values, c=errors, edgecolors="k", s=50, alpha=1, cmap="viridis")
    plt.colorbar(sc, label=r"$L^2$ Error: $\frac{\| u - u_{\mathrm{ref}} \|_2}{\| u_{\mathrm{ref}} \|_2}$")
    plt.xlabel("Mach")
    plt.ylabel("Alpha")
    plt.title("$L^2$ error heat map")
    plt.grid(True, linestyle="--", alpha=0.5)
    # plt.show()
    plt.savefig(f'{plotname}.png', dpi=400)
    plt.close('all')


def plot_error_distribution(datasets={}, plotname='', title='Error distribution'):
    plt.rcParams.update({
    'font.size'            : 16,  # Tamaño general del texto
    'axes.titlesize'       : 20,  # Tamaño del título del gráfico
    'axes.labelsize'       : 18,  # Tamaño de etiquetas de ejes
    'xtick.labelsize'      : 10,  # Tamaño de etiquetas del eje x
    'ytick.labelsize'      : 16,  # Tamaño de etiquetas del eje y
    'legend.fontsize'      : 14,  # Tamaño de la leyenda
    'legend.title_fontsize': 16   # Tamaño del título de la leyenda
    })
    data = []
    dataset_order = []
    dataset_groups = {}
    group_colors = {
        "Train": "lightgreen",
        "Test": "lightgray",
        "Validation   ": "lightblue"
    }
    for dataset_name, (errors, group) in datasets.items():
        data.extend([(dataset_name, err) for err in errors])
        dataset_order.append(dataset_name)
        dataset_groups[dataset_name] = group
    df = pd.DataFrame(data, columns=["Dataset", "Error"])
    plt.figure(figsize=(11, 6))
    # plt.figure()
    ax = plt.gca()  
    for i, dataset_name in enumerate(dataset_order):
        group = dataset_groups[dataset_name]
        color = group_colors.get(group, "white")  
        ax.axvspan(i - 0.5, i + 0.5, facecolor=color, alpha=0.5)
    sns.boxplot(x="Dataset", y="Error", data=df, hue="Dataset", palette="Set2", width=0.1, legend='full',
                showmeans=True, meanprops={"marker": "D", "markerfacecolor": "red", "markeredgecolor": "black", "alpha": 0.5})
    handles, labels = ax.get_legend_handles_labels()
    legend1 = ax.legend(handles, labels, loc='upper left', title="Box Legend", bbox_to_anchor=(1.0,1.02))
    legend_patches = [Patch(facecolor=color, edgecolor='black', alpha=0.5, label=group) 
                      for group, color in group_colors.items()]
    legend2 = plt.legend(handles=legend_patches, loc='upper left', title="Background", bbox_to_anchor=(1.0,0.28))
    ax.add_artist(legend1) 
    ax.xaxis.set_major_locator(MultipleLocator(1))  
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))  
    ax.yaxis.set_major_locator(MultipleLocator(0.05))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))  
    plt.title(f"{title}")
    plt.ylabel(r"Relative $L^2$ Error: $\frac{\| u - u_{\mathrm{ref}} \|_2}{\| u_{\mathrm{ref}} \|_2}$")
    plt.xticks(rotation=0)
    plt.grid(True, linestyle="--", alpha=0.5, which='both')
    plt.tight_layout()
    plt.savefig(f'{plotname}.png', dpi=400)
    plt.close('all')


def rbf_rom_prediction(rom=None, mu_list=[], target_directory='', keyword='', launch_fake_simulation=False, fake_simulation_target_directory=''):
    if not len(mu_list):
        print('There are no elements in mu_list.')
    else:
        #### PREDICTION
        #################################################################################################
        interpolated_solutions = [rom.predict([element]).snapshots_matrix.T for element in mu_list]
        for i, solution in enumerate(interpolated_solutions):
            os.makedirs(f"{target_directory}/RBF_{keyword}", exist_ok=True)
            np.save(f"{target_directory}/RBF_{keyword}/{mu_list[i][0]}, {mu_list[i][1]}", solution)
            #To buid vtk, gid, full_cp and skin_cp from interpolated data and global ids for skin cp
            if launch_fake_simulation:
                LaunchFakeSimulation(solution, [mu_list[i][0], mu_list[i][1], fake_simulation_target_directory, keyword])
            if 'full_cp' in target_directory and os.path.exists('output_plots/RBF/potential_variables/ids_skin_full.npy'):
                directory = 'Skin_Data_skin_cp_from_full_cp'
                ids_skin_full = np.load('output_plots/RBF/potential_variables/ids_skin_full.npy').astype(np.int64)
                os.makedirs(f"{target_directory}/RBF_{directory}", exist_ok=True)
                np.save(f"{target_directory}/RBF_{directory}/{mu_list[i][0]}, {mu_list[i][1]}", solution[ids_skin_full])


def build_rbf(target_directory='output_plots', mu_id=0, weight=1, svd_tol=1e-12, alpha=1.0, beta=1.0, keyword='', QoI='', launch_fake_simulation=False):
    os.makedirs(target_directory, exist_ok=True)
    data_base = initialize_database(mu_id=mu_id,weight=weight,svd_tol=svd_tol,alpha=alpha,beta=beta)   
    mu_train      = [[mu[0],mu[1]] for mu in load_mu_parameters(f'mu_train_{mu_id}_{alpha}_{beta}')]  
    mu_test       = [[mu[0],mu[1]] for mu in load_mu_parameters('mu_test')]  
    mu_validation = [[mu[0],mu[1]] for mu in load_mu_parameters('mu_validation')]  
    if 'potential_variables' in keyword:
        snapshots = data_base.get_snapshots_matrix_from_database(load_mu_parameters(f'mu_train_{mu_id}_{alpha}_{beta}'), table_name='FOM')
    else:
        snapshots = build_snapshots_matrix(data_base=data_base, mu_list=mu_train, mu_id=mu_id, weight=weight, svd_tol=svd_tol, alpha=alpha, beta=beta, QoI=QoI)
    if len(mu_train) < 3:
        print('At least 3 case are Necessary for this method.')
    else:
        parameters = np.array(mu_train)
        rom = initialize_rbf(parameters=parameters, snapshots=snapshots)
        #################################################################################################   
        #### PREDICTION OF TRAIN
        rbf_rom_prediction(rom=rom, mu_list=mu_train, target_directory=target_directory, keyword=keyword,
                           launch_fake_simulation=launch_fake_simulation, fake_simulation_target_directory=target_directory) 
        #################################################################################################
        #### PREDICTION OF TEST
        rbf_rom_prediction(rom=rom, mu_list=mu_test, target_directory=target_directory, keyword=keyword,
                           launch_fake_simulation=launch_fake_simulation, fake_simulation_target_directory=target_directory) 
        #################################################################################################
        #### PREDICTION OF VALIDATION
        rbf_rom_prediction(rom=rom, mu_list=mu_validation, target_directory=target_directory, keyword=keyword,
                           launch_fake_simulation=launch_fake_simulation, fake_simulation_target_directory=target_directory) 
        #################################################################################################


# To save de vtk and GiD models, global position of skin values vector, full and skin cp from dataset
#############################################################################################################################################################
def LaunchFakeSimulation(data_set, mu):
    with open('ProjectParametersNaca0012.json','r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    model = KratosMultiphysics.Model()
    parameters_copy = FakeProjectParameters(parameters.Clone(), mu)
    analysis_stage_class = _GetAnalysisStageClass(parameters_copy)
    simulation = FakeSimulation(analysis_stage_class, model, parameters_copy, data_set, mu)
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
    vtk_output = VtkOutputProcess(simulation.model,
                                    parameters_output)
    vtk_output.ExecuteInitialize()
    vtk_output.ExecuteBeforeSolutionLoop()
    vtk_output.ExecuteInitializeSolutionStep()
    vtk_output.PrintOutput()
    vtk_output.ExecuteFinalizeSolutionStep()
    vtk_output.ExecuteFinalize()
def FakeSimulation(cls, global_model, parameters, data_set, mu):
    class CustomSimulation(cls):
        def __init__(self, model,project_parameters):
            super().__init__(model,project_parameters)
        def Initialize(self):
            super().Initialize()
            model_part = self.model["MainModelPart"]
            if 'full_cp' in mu[3]:
                for (node, value) in zip(model_part.Nodes, data_set):
                    node.SetValue(KratosMultiphysics.REACTION_WATER_PRESSURE, value)
                    node.Fix(KratosMultiphysics.REACTION_WATER_PRESSURE)
            else:
                for node in model_part.Nodes:
                    offset = np.where(np.arange(1,model_part.NumberOfNodes()+1, dtype=int) == node.Id)[0][0]*2
                    node.SetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL, data_set[offset])
                    node.SetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL, data_set[offset+1])
        def Run(self):
            self.Initialize()
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()
            self.Finalize()
        def InitializeSolutionStep(self):
            super().InitializeSolutionStep()
        def FinalizeSolutionStep(self):
            super().FinalizeSolutionStep()
            if 'full_cp' not in mu[3]:
                Cp_final_data = []; ids_skin = []
                modelpart = self.model["MainModelPart.Body2D_Body"]
                for node in modelpart.Nodes:
                    ids_skin.append(node.Id)
                    Cp_final_data.append([node.X, node.Y, node.Z, node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)])
                Full_Cp_final_data = []; ids_skin_full = []
                modelpart = self.model["MainModelPart"]
                for id in ids_skin:
                    for pos, node in enumerate(modelpart.Nodes):
                        if node.Id == id:
                            ids_skin_full.append(pos)
                for node in modelpart.Nodes:
                    Full_Cp_final_data.append([node.X, node.Y, node.Z, node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)])
                os.makedirs(f'{mu[2]}/RBF_Skin_Data_skin_cp_from_potential_variables', exist_ok=True)
                os.makedirs(f'{mu[2]}/RBF_Snapshots_full_cp_from_potential_variables', exist_ok=True)
                np.save(f'{mu[2]}/RBF_Skin_Data_skin_cp_from_potential_variables/{mu[0]}, {mu[1]}', np.array(Cp_final_data))
                np.save(f'{mu[2]}/RBF_Snapshots_full_cp_from_potential_variables/{mu[0]}, {mu[1]}', np.array(Full_Cp_final_data))
                np.save(f'{mu[2]}/ids_skin_full', np.array(ids_skin_full))
    return CustomSimulation(global_model, parameters)
def FakeProjectParameters(parameters, mu):
    angle_of_attack = mu[0]
    mach_infinity   = mu[1]
    parameters["solver_settings"]["echo_level"].SetInt(0)
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(angle_of_attack)
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach_infinity)
    parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString(f'{mu[2]}/Results/{angle_of_attack}, {mach_infinity}')
    parameters["output_processes"]["vtk_output"][0]["Parameters"]["output_path"].SetString(f'{mu[2]}/Results/vtk_output_{angle_of_attack}, {mach_infinity}')
    return parameters
def _GetAnalysisStageClass(parameters):
    analysis_stage_module_name = parameters["analysis_stage"].GetString()
    analysis_stage_class_name = analysis_stage_module_name.split('.')[-1]
    analysis_stage_class_name = ''.join(x.title() for x in analysis_stage_class_name.split('_'))
    analysis_stage_module = importlib.import_module(analysis_stage_module_name)
    analysis_stage_class = getattr(analysis_stage_module, analysis_stage_class_name)
    return analysis_stage_class
#############################################################################################################################################################

def plot_Cps(mu, capture_directory, fom_skin_cp, rbf_skin_cp_from_variables, rbf_skin_cp_from_full_cp, rbf_skin_cp):
    os.makedirs(capture_directory, exist_ok=True)
    case_name = f"{mu[0]}, {mu[1]}"
    capture_filename = f"{capture_directory}/{case_name}.png"
    plt.figure(figsize=(9, 6))
    cp_validation = 0
    # VALIDATION
    if 'validation' in capture_directory:
        validation_filename = f"reference_data/flo36/{case_name}.dat"
        if os.path.exists(validation_filename):
            x = np.loadtxt(validation_filename, usecols=(0,))
            cp_validation = np.loadtxt(validation_filename, usecols=(1,))
            plt.plot(x, cp_validation, '*-', markersize=10.0, label='flo36')
    # FOM
    x_fom, y_fom, cp_fom = fom_skin_cp[:, 0].reshape(-1,1), fom_skin_cp[:, 1].reshape(-1,1), fom_skin_cp[:, 3].reshape(-1,1)
    if 'validation' in capture_directory:
        mask_gtz = y_fom > 0
        mask_lez = ~mask_gtz
        x_sorted_gtz = x_fom[mask_gtz][np.argsort(x_fom[mask_gtz])]
        cp_sorted_gtz = cp_fom[mask_gtz][np.argsort(x_fom[mask_gtz])]
        x_sorted_lez = x_fom[mask_lez][np.argsort(x_fom[mask_lez])]
        cp_sorted_lez = cp_fom[mask_lez][np.argsort(x_fom[mask_lez])]
        cp_val_gtz = np.interp(x_sorted_gtz, x[79:][::-1], cp_validation[79:][::-1])
        cp_val_lez = np.interp(x_sorted_lez, x[:78], cp_validation[:78])
        cp_sim = np.concatenate((cp_sorted_lez, cp_sorted_gtz[::-1]))
        cp_val = np.concatenate((cp_val_lez, cp_val_gtz[::-1]))
        error = np.linalg.norm(cp_sim - cp_val) / np.linalg.norm(cp_val)
        plt.plot(x_fom, cp_fom, 's', markersize=3.5, label=f'FOM-Validation e: {error:.2E}')
    else:
        plt.plot(x_fom, cp_fom, 's', markersize=3.5, label='FOM')
    # RBF Comparisons
    x_rbf, cp_rbf_from_variables = rbf_skin_cp_from_variables[:, 0].reshape(-1,1), rbf_skin_cp_from_variables[:, 3].reshape(-1,1)
    plt.plot(x_rbf, cp_rbf_from_variables, '.', markersize=8.0, label=f'RBF_pot_vars-FOM e: {(np.linalg.norm(cp_fom - cp_rbf_from_variables) / np.linalg.norm(cp_fom)):.2E}')
    plt.plot(x_rbf, rbf_skin_cp_from_full_cp, '.', markersize=6.0, label=f'RBF_full_cp-FOM e: {(np.linalg.norm(cp_fom - rbf_skin_cp_from_full_cp) / np.linalg.norm(cp_fom)):.2E}')
    plt.plot(x_rbf, rbf_skin_cp, '.', markersize=2.0, label=f'RBF_skin_cp-FOM e: {(np.linalg.norm(cp_fom - rbf_skin_cp) / np.linalg.norm(cp_fom)):.2E}')
    cp_min = np.nanmin([np.min(cp_fom), np.min(cp_rbf_from_variables), np.min(rbf_skin_cp_from_full_cp), np.min(rbf_skin_cp), np.min(cp_validation)])
    cp_max = np.nanmax([np.max(cp_fom), np.max(cp_rbf_from_variables), np.max(rbf_skin_cp_from_full_cp), np.max(rbf_skin_cp), np.max(cp_validation)])
    plt.axis([-0.05,1.05,cp_max+0.1,cp_min-0.1])
    plt.title('Cp vs x')
    plt.xlabel('x')
    plt.ylabel('-Cp')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(capture_filename, dpi=200)
    plt.close('all')


def plot_Cps_2d(mu, capture_directory, case_types=['FOM'], rbf_path='output_plots/RBF/potential_variables'):
    os.makedirs(capture_directory, exist_ok=True)
    case_name = f"{mu[0]}, {mu[1]}.png"
    simulation_type = next((s for s in ['Fit', 'Test', 'Run'] if s in capture_directory), 'Unknown')
    for case in case_types:
        path = f"{rbf_path}/Results/vtk_output_{mu[0]}, {mu[1]}" if 'RBF' in case else f"Results/vtk_output_{case}_{simulation_type}{', '.join(map(str, mu))}"
        if os.path.exists(path):
            vtk_filename = get_latest_vtk_file(path)
            if os.path.exists(vtk_filename):
                mesh = pv.read(vtk_filename).extract_geometry()
                points, cells = mesh.points, mesh.faces.reshape(-1, 4)[:, 1:4]
                if 'full_cp' in rbf_path:
                    pressure_coeff = mesh.point_data['REACTION_WATER_PRESSURE']
                else:
                    pressure_coeff = mesh.point_data['PRESSURE_COEFFICIENT']
                fig, ax = plt.subplots(figsize=(10, 8))
                triang = tri.Triangulation(points[:, 0], points[:, 1], cells)
                cmap = cm.get_cmap('viridis') # viridis coolwarm
                contour = ax.tricontour(triang, pressure_coeff, levels=15, cmap='coolwarm', linewidths=0.8)
                ax.tricontourf(triang, pressure_coeff, levels=100, cmap=cmap, alpha=1.0)
                cbar = plt.colorbar(cm.ScalarMappable(norm=plt.Normalize(vmin=pressure_coeff.min(), vmax=pressure_coeff.max()), cmap=cmap), ax=ax, label='PRESSURE COEFFICIENT')
                cbar.ax.tick_params(labelsize=10)
                ax.set_xlabel('X')
                ax.set_ylabel('Y')
                ax.set_title(f'{case} Pressure Coefficient Distribution - Alpha={mu[0]:.3f}, Mach={mu[1]:.3f}', size=13, pad=10)
                ax.set_xlim(-0.5, 1.5)
                ax.set_ylim(-0.75, 1.0)
                plt.savefig(f'{capture_directory}/2D_{case}_{case_name}', dpi=200)
                plt.close('all')


def get_latest_vtk_file(directory):
    vtk_files = [(int(m.group(1)) if m.group(1) else -1, f) for f in os.listdir(directory) if (m := re.match(r"MainModelPart_0(?:_(\d+))?\.vtk", f))]
    return os.path.join(directory, max(vtk_files, key=lambda x: x[0])[1]) if vtk_files else FileNotFoundError("No matching VTK files found.")


def plot_Cps_distributions(mu_list, target_directory, case, **params):
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


def build_snapshots_matrix(data_base=None, mu_list=[], mu_id=0, weight=1, svd_tol=1e-12, alpha=1.0, beta=1.0, QoI='', case='train'):
    snapshots = []
    for mu in mu_list:
        # Necessary to search in the database
        mu = [mu[0],mu[1],weight,beta,mu_id,svd_tol]
        in_database, hash_cp_data = data_base.check_if_in_database("QoI_FOM", mu)
        if not in_database:
            print(f'There is no cp data for mu_{case}: {mu[0]}, {mu[1]}')
            print(f'id={mu_id}, Alpha={alpha}, Beta={beta}.')
        else:
            data = data_base.get_single_numpy_from_database(f'{hash_cp_data}{QoI}')
            snapshots.append(data[:,3].reshape(-1,1))
    return np.block(snapshots)


def build_snapshots_matrix_rbf(mu_list=[], target_directory=''):
    snapshots = []
    for mu in mu_list:
        if np.load(f'{target_directory}/{mu[0]}, {mu[1]}.npy').shape[1] > 1:
            snapshots.append(np.load(f'{target_directory}/{mu[0]}, {mu[1]}.npy')[:,3].reshape(-1,1))
        else:
            snapshots.append(np.load(f'{target_directory}/{mu[0]}, {mu[1]}.npy').reshape(-1,1))
    return np.block(snapshots)


def compute_errors(target_directory='', data_base=None, mu_list=[], mu_id=0, weight=1.0, svd_tol=1e-12, alpha=1.0, beta=1.0, case='', keyword=''):
    if case == 'train':
        mu_parameters = f'mu_train_{mu_id}_{alpha}_{beta}'
    else:
        mu_parameters = f'mu_{case}'
    directories = []; bases = []; QoIs = []
    if 'Snapshots_potential_variables' in keyword:
        directories = [keyword, 'Snapshots_full_cp_from_potential_variables', 'Skin_Data_skin_cp_from_potential_variables']
        bases = ['potential_variables', 'full_cp_from_potential_variables', 'skin_cp_from_potential_variables']
        QoIs = ['', 'Full_Cp_data', 'Cp_data']
    if 'Snapshots_full_cp' in keyword:
        directories = [keyword, 'Skin_Data_skin_cp_from_full_cp']
        bases = ['full_cp', 'skin_cp_from_full_cp']
        QoIs = ['Full_Cp_data', 'Cp_data']
    if 'Skin_Data_skin_cp' in keyword:
        directories = [keyword]   
        bases = ['skin_cp']
        QoIs = ['Cp_data']
    for directory, base, QoI in zip(directories, bases, QoIs):
        if 'cp' not in directory:
            fom_snapshots = data_base.get_snapshots_matrix_from_database(load_mu_parameters(f'{mu_parameters}'), table_name='FOM')
        else:
            fom_snapshots = build_snapshots_matrix(data_base=data_base, mu_list=mu_list, 
                                                mu_id=mu_id, weight=weight, svd_tol=svd_tol, alpha=alpha, beta=beta,
                                                QoI=QoI, case=case)
        rbf_snapshots = build_snapshots_matrix_rbf(mu_list=mu_list, target_directory=f'{target_directory}/RBF_{directory}')
        L2_error = np.linalg.norm(fom_snapshots - rbf_snapshots)/np.linalg.norm(fom_snapshots)
        print('::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::') 
        print(f'RBF {base} approximation error {case}: {L2_error:.2E}')     
        print('::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::') 
        L2_errors = np.linalg.norm(fom_snapshots - rbf_snapshots, axis=0)/np.linalg.norm(fom_snapshots, axis=0)
        np.save(f'{target_directory}/L2_errors_{base}_{case}', L2_errors)


def print_errors(target_directory='output_plots', mu_id=0, weight=1, svd_tol=1e-12, alpha=1.0, beta=1.0, keyword=''):
    data_base = initialize_database(mu_id=mu_id,weight=weight,svd_tol=svd_tol,alpha=alpha,beta=beta) 
    mu_train      = [[mu[0],mu[1]] for mu in load_mu_parameters(f'mu_train_{mu_id}_{alpha}_{beta}')]  
    mu_test       = [[mu[0],mu[1]] for mu in load_mu_parameters('mu_test')]   
    mu_validation = [[mu[0],mu[1]] for mu in load_mu_parameters('mu_validation')] 
    # COMPUTE TRAINING ERRORS
    compute_errors(target_directory=target_directory, data_base=data_base,
                   mu_list=mu_train, mu_id=mu_id, weight=weight, svd_tol=svd_tol, alpha=alpha, beta=beta, case='train', keyword=keyword) 
    # COMPUTE TEST ERRORS
    compute_errors(target_directory=target_directory, data_base=data_base,
                   mu_list=mu_test, mu_id=mu_id, weight=weight, svd_tol=svd_tol, alpha=alpha, beta=beta, case='test', keyword=keyword) 
    # COMPUTE VALIDATION ERRORS
    compute_errors(target_directory=target_directory, data_base=data_base,
                   mu_list=mu_validation, mu_id=mu_id, weight=weight, svd_tol=svd_tol, alpha=alpha, beta=beta, case='validation', keyword=keyword) 




if __name__ == "__main__":
    #Disable logs
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    os.makedirs('output_plots', exist_ok=True)

    # Necessary to define the case to use for the database
    mu_id            = 0      
    weight           = 1.0    
    svd_tol          = 1e-12  
    alpha            = 1.0   
    beta             = 1.0   
    solver_strategy  = 1     
    target_directory = 'output_plots/RBF'
    build_data               = True
    plot_cps_captures        = True
    print_errors_by_models   = True
    plot_heat_maps           = True
    plot_error_distributions = True

    if build_data:
        KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('output_plots/RBF')
    os.makedirs(target_directory, exist_ok=True)

    if build_data:
        #RBF model using all potential flow variables
        build_rbf(target_directory=f'{target_directory}/potential_variables',
                  mu_id=mu_id,weight=weight,svd_tol=svd_tol,alpha=alpha,beta=beta,
                  keyword='Snapshots_potential_variables', launch_fake_simulation=True)
        
        #RBF model using all cp distribution
        build_rbf(target_directory=f'{target_directory}/full_cp',
                    mu_id=mu_id,weight=weight,svd_tol=svd_tol,alpha=alpha,beta=beta,
                    keyword='Snapshots_full_cp',QoI='Full_Cp_data', launch_fake_simulation=True)
        
        # RBF model using only cp's of the skin
        build_rbf(target_directory=f'{target_directory}/skin_cp',
                    mu_id=mu_id,weight=weight,svd_tol=svd_tol,alpha=alpha,beta=beta,
                    keyword='Skin_Data_skin_cp',QoI='Cp_data')

    # Posprocess -- Plot Cp distributions
    if plot_cps_captures:
        directories = ['Training', 'Test', 'Validation']
        for directory in directories:
            KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting(f'{target_directory}/{directory}')
            os.makedirs(f'{target_directory}/{directory}', exist_ok=True)

        mu_train      = [[mu[0],mu[1]] for mu in load_mu_parameters(f'mu_train_{mu_id}_{alpha}_{beta}')]  
        mu_test       = [[mu[0],mu[1]] for mu in load_mu_parameters('mu_test')]   
        mu_validation = [[mu[0],mu[1]] for mu in load_mu_parameters('mu_validation')]  

        plot_Cps_distributions(mu_list=mu_train, target_directory=target_directory,
                            case='Fit', mu_id=mu_id, weight=weight, svd_tol=svd_tol, alpha=alpha, beta=beta, solver_strategy=solver_strategy)

        plot_Cps_distributions(mu_list=mu_test, target_directory=target_directory,
                            case='Test', mu_id=mu_id, weight=weight, svd_tol=svd_tol, alpha=alpha, beta=beta, solver_strategy=solver_strategy)

        plot_Cps_distributions(mu_list=mu_validation, target_directory=target_directory,
                            case='Run_validation', mu_id=mu_id, weight=weight, svd_tol=svd_tol, alpha=alpha, beta=beta, solver_strategy=solver_strategy)

    # Posprocess -- Print errors
    if print_errors_by_models:
        print_errors(target_directory=f'{target_directory}/potential_variables',
                    mu_id=mu_id, weight=weight, svd_tol=svd_tol, alpha=alpha, beta=beta,
                    keyword='Snapshots_potential_variables')
        
        print_errors(target_directory=f'{target_directory}/full_cp',
                    mu_id=mu_id, weight=weight, svd_tol=svd_tol, alpha=alpha, beta=beta,
                    keyword='Snapshots_full_cp')
        
        print_errors(target_directory=f'{target_directory}/skin_cp',
                    mu_id=mu_id, weight=weight, svd_tol=svd_tol, alpha=alpha, beta=beta,
                    keyword='Skin_Data_skin_cp')

    # Posprocess -- Plot L2 error heatmaps
    if plot_heat_maps:
        mu_train      = [[mu[0],mu[1]] for mu in load_mu_parameters(f'mu_train_{mu_id}_{alpha}_{beta}')]  
        mu_test       = [[mu[0],mu[1]] for mu in load_mu_parameters('mu_test')]   
        mu_validation = [[mu[0],mu[1]] for mu in load_mu_parameters('mu_validation')] 

        # LOAD ERRORS
        # POTENTIAL VARIABLES
        #################################################################################################################################
        # TRAINING
        error_var_tr       = np.load(f'{target_directory}/potential_variables/L2_errors_potential_variables_train.npy')
        error_fcp_f_var_tr = np.load(f'{target_directory}/potential_variables/L2_errors_full_cp_from_potential_variables_train.npy')
        error_scp_f_var_tr = np.load(f'{target_directory}/potential_variables/L2_errors_skin_cp_from_potential_variables_train.npy')
        # TEST
        error_var_t        = np.load(f'{target_directory}/potential_variables/L2_errors_potential_variables_test.npy')
        error_fcp_f_var_t  = np.load(f'{target_directory}/potential_variables/L2_errors_full_cp_from_potential_variables_test.npy')
        error_scp_f_var_t  = np.load(f'{target_directory}/potential_variables/L2_errors_skin_cp_from_potential_variables_test.npy')
        # VALIDATION
        error_var_v        = np.load(f'{target_directory}/potential_variables/L2_errors_potential_variables_validation.npy')
        error_fcp_f_var_v  = np.load(f'{target_directory}/potential_variables/L2_errors_full_cp_from_potential_variables_validation.npy')
        error_scp_f_var_v  = np.load(f'{target_directory}/potential_variables/L2_errors_skin_cp_from_potential_variables_validation.npy')
        # FULL CP
        #################################################################################################################################
        # TRAINING
        error_fcp_tr       = np.load(f'{target_directory}/full_cp/L2_errors_full_cp_train.npy')
        error_scp_f_fcp_tr = np.load(f'{target_directory}/full_cp/L2_errors_skin_cp_from_full_cp_train.npy')
        # TEST
        error_fcp_t        = np.load(f'{target_directory}/full_cp/L2_errors_full_cp_test.npy')
        error_scp_f_fcp_t  = np.load(f'{target_directory}/full_cp/L2_errors_skin_cp_from_full_cp_test.npy')
        # VALIDATION
        error_fcp_v        = np.load(f'{target_directory}/full_cp/L2_errors_full_cp_validation.npy')
        error_scp_f_fcp_v  = np.load(f'{target_directory}/full_cp/L2_errors_skin_cp_from_full_cp_validation.npy')
        # SKIN CP
        #################################################################################################################################
        # TRAINING
        error_scp_tr       = np.load(f'{target_directory}/skin_cp/L2_errors_skin_cp_train.npy')
        # TEST
        error_scp_t        = np.load(f'{target_directory}/skin_cp/L2_errors_skin_cp_test.npy')
        # VALIDATION
        error_scp_v        = np.load(f'{target_directory}/skin_cp/L2_errors_skin_cp_validation.npy')
        #################################################################################################################################

        plot_L2_error_heat_map_interpolated(plotname=f'{target_directory}/L2_error_skin_cp_from_variables_heat_map_interpolated_test', 
                                            mu_list=mu_test, errors=error_scp_f_var_t)
        
        plot_L2_error_heat_map_discreted(plotname=f'{target_directory}/L2_error_skin_cp_from_variables_heat_map_discretet_test', 
                                            mu_list=mu_test, errors=error_scp_f_var_t)

    # Posprocess -- Plot L2 error distributions
    if plot_error_distributions:
        # POTENTIAL VARIABLES
        #################################################################################################################################
        # TRAINING
        error_var_tr       = np.load(f'{target_directory}/potential_variables/L2_errors_potential_variables_train.npy')
        error_fcp_f_var_tr = np.load(f'{target_directory}/potential_variables/L2_errors_full_cp_from_potential_variables_train.npy')
        error_scp_f_var_tr = np.load(f'{target_directory}/potential_variables/L2_errors_skin_cp_from_potential_variables_train.npy')
        # TEST
        error_var_t        = np.load(f'{target_directory}/potential_variables/L2_errors_potential_variables_test.npy')
        error_fcp_f_var_t  = np.load(f'{target_directory}/potential_variables/L2_errors_full_cp_from_potential_variables_test.npy')
        error_scp_f_var_t  = np.load(f'{target_directory}/potential_variables/L2_errors_skin_cp_from_potential_variables_test.npy')
        # VALIDATION
        error_var_v        = np.load(f'{target_directory}/potential_variables/L2_errors_potential_variables_validation.npy')
        error_fcp_f_var_v  = np.load(f'{target_directory}/potential_variables/L2_errors_full_cp_from_potential_variables_validation.npy')
        error_scp_f_var_v  = np.load(f'{target_directory}/potential_variables/L2_errors_skin_cp_from_potential_variables_validation.npy')
        # FULL CP
        #################################################################################################################################
        # TRAINING
        error_fcp_tr       = np.load(f'{target_directory}/full_cp/L2_errors_full_cp_train.npy')
        error_scp_f_fcp_tr = np.load(f'{target_directory}/full_cp/L2_errors_skin_cp_from_full_cp_train.npy')
        # TEST
        error_fcp_t        = np.load(f'{target_directory}/full_cp/L2_errors_full_cp_test.npy')
        error_scp_f_fcp_t  = np.load(f'{target_directory}/full_cp/L2_errors_skin_cp_from_full_cp_test.npy')
        # VALIDATION
        error_fcp_v        = np.load(f'{target_directory}/full_cp/L2_errors_full_cp_validation.npy')
        error_scp_f_fcp_v  = np.load(f'{target_directory}/full_cp/L2_errors_skin_cp_from_full_cp_validation.npy')
        # SKIN CP
        #################################################################################################################################
        # TRAINING
        error_scp_tr       = np.load(f'{target_directory}/skin_cp/L2_errors_skin_cp_train.npy')
        # TEST
        error_scp_t        = np.load(f'{target_directory}/skin_cp/L2_errors_skin_cp_test.npy')
        # VALIDATION
        error_scp_v        = np.load(f'{target_directory}/skin_cp/L2_errors_skin_cp_validation.npy')
        #################################################################################################################################

        datasets = {
            "Potential_tr": (error_scp_f_var_tr, 'Train'),
            "Full Cp_tr"  : (error_scp_f_fcp_tr, 'Train'),
            "Skin Cp_tr"  : (error_scp_tr, 'Train'),
            "Potential_t ": (error_scp_f_var_t, 'Test'),
            "Full Cp_t "  : (error_scp_f_fcp_t, 'Test'),
            "Skin Cp_t "  : (error_scp_t, 'Test'),
            "Potential_v ": (error_scp_f_var_v, 'Validation   '),
            "Full Cp_v "  : (error_scp_f_fcp_v, 'Validation   '),
            "Skin Cp_v "  : (error_scp_v, 'Validation   ')
        }

        plot_error_distribution(datasets=datasets, 
                                plotname=f'{target_directory}/Error_distribution_skin_data', 
                                title='Error distribution - skin Cp')