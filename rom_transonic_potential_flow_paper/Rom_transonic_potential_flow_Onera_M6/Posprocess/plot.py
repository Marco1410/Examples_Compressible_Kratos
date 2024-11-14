import os
import importlib
import numpy as np
import pyvista as pv 
import KratosMultiphysics
import matplotlib.cm as cm
import matplotlib.tri as tri
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from KratosMultiphysics.gid_output_process import GiDOutputProcess
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp


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
            modelpart = self.model["MainModelPart.Wing"]
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
    parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString(f'Results/{mu[2]}_{angle_of_attack}, {mach_infinity}')
    parameters["output_processes"]["vtk_output"][0]["Parameters"]["output_path"].SetString(f'Results/{mu[2]}_{angle_of_attack}, {mach_infinity}')
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
def load_mu_parameters(name):
    filename = f'mu_{name}.npy'
    if os.path.exists(filename):
        mu_npy = np.load(filename)
        mu =  [mu.tolist() for mu in mu_npy]
    else:
        mu = []
    return mu

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# save parameters
#
def Plot_Cps(mu_list, capture_directory, plot=['CP']):
    #Disable logs
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    for id, mu in enumerate(mu_list):
        case_name = f'{mu[0]}, {mu[1]}'
        capture_filename   = f"{case_name}.png"

        ########## 3D PLOT #####################################################################
        if '3D' in plot:
            vtk_filename = f"Results/vtk_output_ROM_Fit{id}/MainModelPart_Wing_0_1.vtk"
            if os.path.exists(vtk_filename):
                ######################################################################
                # Cargar el archivo .vtk
                mesh = pv.read(vtk_filename)

                # Extraer la geometría superficial
                surface_mesh = mesh.extract_geometry()

                # Extraer los puntos y las celdas
                points = surface_mesh.points
                cells = surface_mesh.faces.reshape(-1, 4)[:, 1:4]

                # Extraer los valores de la variable 'PRESSURE_COEFFICIENT'
                pressure_coeff = surface_mesh.point_data['PRESSURE_COEFFICIENT']

                # Normalizar los valores del coeficiente de presión para mapearlos a la escala de colores
                norm = plt.Normalize(vmin=pressure_coeff.min(), vmax=pressure_coeff.max())
                cmap = cm.get_cmap('viridis')

                ax = plt.figure(figsize=(13, 10)).add_subplot(projection='3d')

                triangles = [points[cell] for cell in cells]
                facecolors = []

                for cell in cells:
                    vertex_colors = cmap(norm(pressure_coeff[cell]))
                    facecolors.append(np.mean(vertex_colors, axis=0))

                tri_collection = Poly3DCollection(triangles,
                                                facecolors=facecolors, 
                                                linewidth=0.2, 
                                                edgecolor=facecolors,
                                                antialiased=True, 
                                                alpha=1.0)

                ax.add_collection3d(tri_collection)
                ax.view_init(elev=30, azim=135) 

                ax.set_xlim(points[:, 0].min(), points[:, 0].max())
                ax.set_ylim(points[:, 1].min(), points[:, 1].max())
                ax.set_zlim([-0.25, 0.25])

                mappable = cm.ScalarMappable(norm=norm, cmap='viridis')
                mappable.set_array(pressure_coeff)
                plt.colorbar(mappable, ax=ax, label='PRESSURE COEFFICIENT', location='right')

                ax.set_xlabel('X')
                ax.set_ylabel('Y')
                ax.set_zlabel('Z')
                ax.set_title(f' ROM Pressure distribution - Angle={np.round(mu[0],3)}, Mach={np.round(mu[1],3)}', size=13, pad=8, y=-0.15)
                # plt.show()
                plt.savefig(f'{capture_directory}/3D_{capture_filename}', dpi=200)
                plt.close('all')

        ########## 2D PLOT #####################################################################
        if '2D' in plot:
            vtk_filename = f"Results/vtk_output_ROM_Fit{id}/MainModelPart_Wing_0_1.vtk"
            if os.path.exists(vtk_filename):
                ######################################################################
                # Cargar el archivo .vtk
                mesh = pv.read(vtk_filename)

                # Extraer la geometría superficial
                surface_mesh = mesh.extract_geometry()

                # Extraer los puntos y las celdas
                points = surface_mesh.points
                cells = surface_mesh.faces.reshape(-1, 4)[:, 1:4]

                # Extraer los valores de la variable 'PRESSURE_COEFFICIENT'
                pressure_coeff = surface_mesh.point_data['PRESSURE_COEFFICIENT']

                # Normalizar los valores del coeficiente de presión para mapearlos a la escala de colores
                norm = plt.Normalize(vmin=pressure_coeff.min(), vmax=pressure_coeff.max())
                cmap = cm.get_cmap('viridis')

                ax = plt.figure(figsize=(10, 10)).add_subplot()
                # EXTRACT Y > 0 POINTS
                mask = points[:, 1] > 0
                filtered_points = points[mask]
                filtered_pressure_coeff = pressure_coeff[mask]

                index_map = {old_idx: new_idx for new_idx, old_idx in enumerate(np.where(mask)[0])}

                filtered_cells = []
                for cell in cells:
                    if all(mask[cell]):  
                        filtered_cells.append([index_map[idx] for idx in cell])

                triang = tri.Triangulation(filtered_points[:, 0], filtered_points[:, 1], filtered_cells)
                ax.tripcolor(triang, filtered_pressure_coeff, cmap='viridis', linewidth=0.5)

                mappable = cm.ScalarMappable(norm=norm, cmap='viridis')
                mappable.set_array(pressure_coeff)
                plt.colorbar(mappable, ax=ax, label='PRESSURE COEFFICIENT', location='right')
                ax.set_xlabel('X')
                ax.set_ylabel('Y')
                ax.set_title(f' ROM Pressure distribution - Angle={np.round(mu[0],3)}, Mach={np.round(mu[1],3)}', size=13, pad=8, y=-0.125)
                # plt.show()
                plt.savefig(f'{capture_directory}/2D_{capture_filename}', dpi=200)
                plt.close('all')

        ########## CP PLOT #####################################################################
        if 'CP' in plot:
            fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(13, 7))
            labels = [95, 90, 80, 65, 44, 20] 
            fv_potential_lines_in_data_files = [117, 116, 116, 116, 119, 122]

            # selecciona que modelos mostrar
            #########################################################
            show_all     = True
            #########################################################
            VALIDATION   = True
            EXPERIMENTAL = True
            POTENTIAL    = True
            FOM          = False
            ROM          = False
            HROM         = False
            HHROM        = False
            RBF          = False
            #########################################################
            if show_all:
                VALIDATION=FOM=ROM=HROM=HHROM=RBF=True

            #### FOM ######
            fom_skin_data_filename = f"FOM_Skin_Data/{case_name}.dat"
            if os.path.exists(fom_skin_data_filename) and FOM:
                x  = np.loadtxt(fom_skin_data_filename, usecols=(0,))
                y  = np.loadtxt(fom_skin_data_filename, usecols=(1,))
                z  = np.loadtxt(fom_skin_data_filename, usecols=(2,))
                cp_f = np.loadtxt(fom_skin_data_filename, usecols=(3,))

                indexes    = z > 0
                x_fom      = x[indexes]
                y_fom      = y[indexes]
                cp_fom     = cp_f[indexes]
                indexes    = z <= 0
                x_fom_inf  = x[indexes]
                y_fom_inf  = y[indexes]
                cp_fom_inf = cp_f[indexes]

            #### ROM ######
            rom_skin_data_filename = f"ROM_Skin_Data/{case_name}.dat"
            if os.path.exists(rom_skin_data_filename) and ROM:
                x  = np.loadtxt(rom_skin_data_filename, usecols=(0,))
                y  = np.loadtxt(rom_skin_data_filename, usecols=(1,))
                z  = np.loadtxt(rom_skin_data_filename, usecols=(2,))
                cp_r = np.loadtxt(rom_skin_data_filename, usecols=(3,))

                indexes    = z > 0
                x_rom      = x[indexes]
                y_rom      = y[indexes]
                cp_rom     = cp_r[indexes]
                indexes    = z <= 0
                x_rom_inf  = x[indexes]
                y_rom_inf  = y[indexes]
                cp_rom_inf = cp_r[indexes]

            #### HROM ######
            hrom_skin_data_filename = f"HROM_Skin_Data/{case_name}.dat"
            if os.path.exists(hrom_skin_data_filename) and HROM:
                x  = np.loadtxt(hrom_skin_data_filename, usecols=(0,))
                y  = np.loadtxt(hrom_skin_data_filename, usecols=(1,))
                z  = np.loadtxt(hrom_skin_data_filename, usecols=(2,))
                cp_h = np.loadtxt(hrom_skin_data_filename, usecols=(3,))

                indexes    = z > 0
                x_hrom      = x[indexes]
                y_hrom      = y[indexes]
                cp_hrom     = cp_h[indexes]
                indexes    = z <= 0
                x_hrom_inf  = x[indexes]
                y_hrom_inf  = y[indexes]
                cp_hrom_inf = cp_h[indexes]

            #### HHROM ######
            hhrom_name = f"HHROM_Snapshots/{case_name}.npy"
            hhrom_skin_data_filename = f"HHROM_Skin_Data/{case_name}.dat"
            if os.path.exists(hhrom_name) and HHROM:
                hhrom = np.load(hhrom_name)
                LaunchFakeSimulation(hhrom, [mu[0], mu[1], 'HHROM'])
                x  = np.loadtxt(hhrom_skin_data_filename, usecols=(0,))
                y  = np.loadtxt(hhrom_skin_data_filename, usecols=(1,))
                z  = np.loadtxt(hhrom_skin_data_filename, usecols=(2,))
                cp_hh = np.loadtxt(hhrom_skin_data_filename, usecols=(3,))

                indexes    = z > 0
                x_hhrom      = x[indexes]
                y_hhrom      = y[indexes]
                cp_hhrom     = cp_hh[indexes]
                indexes    = z <= 0
                x_hhrom_inf  = x[indexes]
                y_hhrom_inf  = y[indexes]
                cp_hhrom_inf = cp_hh[indexes]

            #### RBF ######
            rbf_skin_data_filename = f"RBF_Skin_Data/{case_name}.npy"
            if os.path.exists(rbf_skin_data_filename) and RBF:
                x  = np.loadtxt(f"FOM_Skin_Data/{case_name}.dat", usecols=(0,))
                y  = np.loadtxt(f"FOM_Skin_Data/{case_name}.dat", usecols=(1,))
                z  = np.loadtxt(f"FOM_Skin_Data/{case_name}.dat", usecols=(2,))
                cp_rb = np.array(np.load(rbf_skin_data_filename)).reshape(-1,1) 

                indexes    = z > 0
                x_rbf      = x[indexes]
                y_rbf      = y[indexes]
                cp_rbf     = cp_rb[indexes]
                indexes    = z <= 0
                x_rbf_inf  = x[indexes]
                y_rbf_inf  = y[indexes]
                cp_rbf_inf = cp_rb[indexes]
            
            if os.path.exists(fom_skin_data_filename) and os.path.exists(rom_skin_data_filename):
                error_rom_fom = np.linalg.norm(cp_f-cp_r)/np.linalg.norm(cp_f)
            if os.path.exists(fom_skin_data_filename) and os.path.exists(hrom_skin_data_filename):
                error_hrom_fom = np.linalg.norm(cp_f-cp_h)/np.linalg.norm(cp_f)
            if os.path.exists(fom_skin_data_filename) and os.path.exists(hhrom_skin_data_filename):
                error_hhrom_fom = np.linalg.norm(cp_f-cp_hh)/np.linalg.norm(cp_f)
            if os.path.exists(fom_skin_data_filename) and os.path.exists(rbf_skin_data_filename):
                error_rbf_fom = np.linalg.norm(np.array(cp_f).reshape(-1,1)-cp_rb)/np.linalg.norm(np.array(cp_f).reshape(-1,1))

            for idx, sub_ax in enumerate(axs.flat):   
                # SET Y POSITION TO PLOT CP DISTRIBUTION
                b = 1.1963
                y_target = labels[idx]*b/100 # y_target = y/b

                # CREATE X LINEAR SPACE
                x_min = (y_target)*np.tan(30*np.pi/180)
                x_max = (y_target)*np.tan(15.8*np.pi/180)+0.8059
                x_grid = np.linspace(x_min, x_max, 150)
                # SET X BETWEEN 0 and 1
                x_airfoil_normalized = (x_grid - x_min) / (x_max - x_min)
                x_airfoil_normalized_full = np.concatenate((x_airfoil_normalized, x_airfoil_normalized))
                
                #### VALIDATION ######
                if capture_directory == 'Validation' and VALIDATION:
                    experimental_skin_data_filename = f'../reference_data/onera/experiment/cp_{labels[idx]}.dat'
                    if os.path.exists(experimental_skin_data_filename) and EXPERIMENTAL:
                        x_validation_exp  = np.loadtxt(experimental_skin_data_filename, usecols=(0,))
                        cp_validation_exp = np.loadtxt(experimental_skin_data_filename, usecols=(1,))
                        # ORIGINAL
                        # sub_ax.scatter(x_validation_exp, cp_validation_exp, marker="+", label="Experimental orig", s=7)
                        ####################################################################################################################
                        # INTERPOLATED
                        cp_interpolated_exp = np.interp(x_airfoil_normalized, x_validation_exp[13:], cp_validation_exp[13:])
                        cp_interpolated_exp_inf = np.interp(x_airfoil_normalized, x_validation_exp[:12][::-1], cp_validation_exp[:12][::-1])
                        cp_validation_exp_interpolated = np.concatenate((cp_interpolated_exp_inf, cp_interpolated_exp))
                        sub_ax.scatter(x_airfoil_normalized_full, cp_validation_exp_interpolated, marker="+", label="Experimental", s=7)

                    #### POTENTIAL ######
                    potential_skin_data_filename = f'../reference_data/onera/potential_solver/cp_{labels[idx]}_new.dat'
                    if os.path.exists(potential_skin_data_filename) and POTENTIAL:
                        x_validation_pot  = np.loadtxt(potential_skin_data_filename, usecols=(0,))
                        cp_validation_pot = np.loadtxt(potential_skin_data_filename, usecols=(1,))
                        # ORIGINAL
                        # sub_ax.scatter(x_validation_pot, cp_validation_pot, marker="*", label="Validation potential solver orig", s=4)
                        ####################################################################################################################
                        # INTERPOLATED
                        cp_interpolated_pot = np.interp(x_airfoil_normalized, x_validation_pot[fv_potential_lines_in_data_files[idx]+1:], cp_validation_pot[fv_potential_lines_in_data_files[idx]+1:])
                        cp_interpolated_pot_inf = np.interp(x_airfoil_normalized, x_validation_pot[:fv_potential_lines_in_data_files[idx]], cp_validation_pot[:fv_potential_lines_in_data_files[idx]])
                        cp_validation_pot_interpolated = np.concatenate((cp_interpolated_pot_inf, cp_interpolated_pot))
                        sub_ax.scatter(x_airfoil_normalized_full, cp_validation_pot_interpolated, marker="*", label="FV Potential solver", s=4)

                #### FOM ######
                if os.path.exists(fom_skin_data_filename) and FOM:
                    # INTERPOLATED
                    cp_interpolated_fom = griddata((x_fom, y_fom), cp_fom, (x_grid, y_target), method='linear', fill_value=0.25)
                    cp_interpolated_fom_inf = griddata((x_fom_inf, y_fom_inf), cp_fom_inf, (x_grid, y_target), method='linear')
                    cp_interpolated_fom_full = np.concatenate((cp_interpolated_fom_inf, cp_interpolated_fom))
                    if capture_directory == 'Validation' and VALIDATION and POTENTIAL:
                        error_pot_fom = np.linalg.norm(-cp_interpolated_fom_full-cp_validation_pot_interpolated)/np.linalg.norm(cp_validation_pot_interpolated)
                        sub_ax.scatter(x_airfoil_normalized_full, -cp_interpolated_fom_full, marker="s", label=f"FOM-FV solver e: {error_pot_fom:.2E}", s=4)
                    else:
                        sub_ax.scatter(x_airfoil_normalized_full, -cp_interpolated_fom_full, marker="s", label=f"FOM", s=4)

                #### ROM ######
                if os.path.exists(rom_skin_data_filename) and ROM:
                    cp_interpolated_rom = griddata((x_rom, y_rom), cp_rom, (x_grid, y_target), method='linear')
                    cp_interpolated_rom_inf = griddata((x_rom_inf, y_rom_inf), cp_rom_inf, (x_grid, y_target), method='linear')
                    cp_interpolated_rom_full = np.concatenate((cp_interpolated_rom_inf, cp_interpolated_rom))
                    sub_ax.scatter(x_airfoil_normalized_full, -cp_interpolated_rom_full, marker="o", label=f"ROM-FOM e: {error_rom_fom:.2E}", s=3)

                #### HROM ######
                if os.path.exists(hrom_skin_data_filename) and HROM:
                    cp_interpolated_hrom = griddata((x_hrom, y_hrom), cp_hrom, (x_grid, y_target), method='linear')
                    cp_interpolated_hrom_inf = griddata((x_hrom_inf, y_hrom_inf), cp_hrom_inf, (x_grid, y_target), method='linear')
                    cp_interpolated_hrom_full = np.concatenate((cp_interpolated_hrom_inf, cp_interpolated_hrom))
                    sub_ax.scatter(x_airfoil_normalized_full, -cp_interpolated_hrom_full, marker="d", label=f"HROM-FOM e: {error_hrom_fom:.2E}", s=6)

                #### HHROM ######
                if os.path.exists(hhrom_name) and HHROM:
                    cp_interpolated_hhrom = griddata((x_hhrom, y_hhrom), cp_hhrom, (x_grid, y_target), method='linear')
                    cp_interpolated_hhrom_inf = griddata((x_hhrom_inf, y_hhrom_inf), cp_hhrom_inf, (x_grid, y_target), method='linear')
                    cp_interpolated_hhrom_full = np.concatenate((cp_interpolated_hhrom_inf, cp_interpolated_hhrom))
                    sub_ax.scatter(x_airfoil_normalized_full, -cp_interpolated_hhrom_full, marker="X", label=f"HHROM-FOM e: {error_hhrom_fom:.2E}", s=6)

                #### RBF ######
                if os.path.exists(rbf_skin_data_filename) and RBF:
                    cp_interpolated_rbf = griddata((x_rbf, y_rbf), cp_rbf, (x_grid, y_target), method='linear')
                    cp_interpolated_rbf_inf = griddata((x_rbf_inf, y_rbf_inf), cp_rbf_inf, (x_grid, y_target), method='linear')
                    cp_interpolated_rbf_full = np.concatenate((cp_interpolated_rbf_inf, cp_interpolated_rbf))
                    sub_ax.scatter(x_airfoil_normalized_full, -cp_interpolated_rbf_full, marker="P", label=f"RBF-FOM e: {error_rbf_fom:.2E}", s=2)

                sub_ax.set_title(f'Section y/b = {np.round(y_target/b,2)}')
                sub_ax.set_xlabel('x')
                sub_ax.set_ylabel('-Cp')

                if idx == 2:
                    sub_ax.legend(bbox_to_anchor=(1.04, 1), loc='upper left', borderaxespad=0.0)
                sub_ax.set_title(f'Section y/b = {np.round(y_target/b,2)}')
                if idx == 3 or idx == 4 or idx == 5:
                    sub_ax.set_xlabel('x')
                if idx == 0 or idx == 3:
                    sub_ax.set_ylabel('-Cp')

            plt.tight_layout()
            # plt.show()
            plt.savefig(f'{capture_directory}/2D_CP_{capture_filename}', dpi=200)
            plt.close('all')


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


if __name__ == "__main__":

    mu_train      = load_mu_parameters('train')
    mu_test       = load_mu_parameters('test')
    mu_validation = load_mu_parameters('validation')

    Plot_Cps(mu_train     , 'Train_Captures', plot=['2D'])
    Plot_Cps(mu_test      , 'Test_Captures' , plot=['2D'])
    Plot_Cps(mu_validation, 'Validation'    , plot=['2D'])
