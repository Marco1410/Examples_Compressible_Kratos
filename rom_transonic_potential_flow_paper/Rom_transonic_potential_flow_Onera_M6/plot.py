import pyvista as pv 
import numpy as np
import os
import importlib
import KratosMultiphysics
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.cm as cm
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
def Plot_Cps(mu_list, capture_directory):
    #Disable logs
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

    for id, mu in enumerate(mu_list):
        case_name = f'{mu[0]}, {mu[1]}'
        capture_filename   = f"{capture_directory}/{case_name}.png"

        #### CP PLOT
        ######################################################################
        # Cargar el archivo .vtk
        mesh = pv.read(f"Results/vtk_output_ROM_Fit{id}/MainModelPart_Wing_0_1.vtk")

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

        # Crear la figura y el eje 3D
        fig = plt.figure(figsize=(12, 10))
        # Posición y tamaño del gráfico principal (ajustable)
        principal_plot_position = [0.05, 0.08, 0.6, 0.6]  # [left, bottom, width, height]
        ax = fig.add_axes(principal_plot_position, projection='3d')

        # Crear los triángulos y los colores correspondientes
        triangles = [points[cell] for cell in cells]
        facecolors = []

        for cell in cells:
            # Obtener los colores para cada vértice del triángulo
            vertex_colors = cmap(norm(pressure_coeff[cell]))
            # Crear un color promedio para el triángulo
            facecolors.append(np.mean(vertex_colors, axis=0))

        # Crear la colección de polígonos con colores interpolados
        tri_collection = Poly3DCollection(triangles,
                                        facecolors=facecolors, 
                                        linewidth=0, 
                                        edgecolor=facecolors,
                                        antialiased=False, 
                                        alpha=1.0)

        # Añadir la colección al gráfico
        ax.add_collection3d(tri_collection)
        ax.view_init(elev=30, azim=135) 

        # Configurar los límites del gráfico para que se ajusten a la malla
        ax.set_xlim(points[:, 0].min(), points[:, 0].max())
        ax.set_ylim(points[:, 1].min(), points[:, 1].max())
        ax.set_zlim([-0.25, 0.25])

        # Añadir una barra de color para la escala del coeficiente de presión
        mappable = cm.ScalarMappable(norm=norm, cmap='viridis')
        mappable.set_array(pressure_coeff)
        plt.colorbar(mappable, ax=ax, label='PRESSURE COEFFICIENT', location='left')

        # Configurar etiquetas y título
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title(f' ROM Pressure distribution - Angle={np.round(mu[0],3)}, Mach={np.round(mu[1],3)}', size=13, pad=8, y=-0.15)

        # Dibujar subplots
        #########################################################
        # Posiciones relativas de los subplots
        positions = [
            [0.08, 0.78, 0.18, 0.18],
            [0.34, 0.78, 0.18, 0.18],
            [0.60, 0.78, 0.18, 0.18],
            [0.80, 0.57, 0.18, 0.18],
            [0.80, 0.31, 0.18, 0.18],
            [0.80, 0.055, 0.18, 0.18]
        ]
        labels = [95, 90, 80, 65, 44, 20] 
        potential_lines = [117, 116, 116, 116, 119, 122]

        # selecciona que modelos mostrar
        show_all     = True
        VALIDATION   = True
        EXPERIMENTAL = True
        POTENTIAL    = True
        FOM          = False
        ROM          = False
        HROM         = False
        HHROM        = False
        RBF          = False
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
            error_rom_fom = np.linalg.norm(cp_f-cp_r)/np.linalg.norm(cp_f)

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
            error_hrom_fom = np.linalg.norm(cp_f-cp_h)/np.linalg.norm(cp_f)

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
        if os.path.exists(hhrom_name) and HHROM:
            hhrom = np.load(hhrom_name)
            LaunchFakeSimulation(hhrom, [mu[0], mu[1], 'HHROM'])
            hhrom_skin_data_filename = f"HHROM_Skin_Data/{case_name}.dat"
            x  = np.loadtxt(hhrom_skin_data_filename, usecols=(0,))
            y  = np.loadtxt(hhrom_skin_data_filename, usecols=(1,))
            z  = np.loadtxt(hhrom_skin_data_filename, usecols=(2,))
            cp_hh = np.loadtxt(hhrom_skin_data_filename, usecols=(3,))
            error_hhrom_fom = np.linalg.norm(cp_f-cp_hh)/np.linalg.norm(cp_f)

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
            error_rbf_fom = np.linalg.norm(np.array(cp_f).reshape(-1,1)-cp_rb)/np.linalg.norm(np.array(cp_f).reshape(-1,1))

            indexes    = z > 0
            x_rbf      = x[indexes]
            y_rbf      = y[indexes]
            cp_rbf     = cp_rb[indexes]
            indexes    = z <= 0
            x_rbf_inf  = x[indexes]
            y_rbf_inf  = y[indexes]
            cp_rbf_inf = cp_rb[indexes]

        # Dibujar subplots
        for idx in range(len(positions)):
            sub_ax = fig.add_axes(positions[idx])    

            # Define la posición específica en el eje Y donde se desea interpolar
            y_target = labels[idx]/100 

            # Crea un espacio lineal de valores X para interpolar
            x_min = (y_target)*np.tan(30*np.pi/180)
            x_max = (y_target)*np.tan(15.8*np.pi/180)+0.8059
            x_grid = np.linspace(x_min, x_max, 150)
            # Normalizar X entre 0 y 1
            x_airfoil_normalized = (x_grid - x_min) / (x_max - x_min)
            x_airfoil_normalized_full = np.concatenate((x_airfoil_normalized, x_airfoil_normalized))
            
            #### VALIDATION ######
            if capture_directory == 'Validation' and VALIDATION:
                experimental_skin_data_filename = f'../reference_data/onera/experiment/cp_{labels[idx]}.dat'
                if os.path.exists(experimental_skin_data_filename) and EXPERIMENTAL:
                    x_validation_exp  = np.loadtxt(experimental_skin_data_filename, usecols=(0,))
                    cp_validation_exp = np.loadtxt(experimental_skin_data_filename, usecols=(1,))
                    cp_interpolated_exp = np.interp(x_airfoil_normalized, x_validation_exp[13:], cp_validation_exp[13:])
                    cp_interpolated_exp_inf = np.interp(x_airfoil_normalized, x_validation_exp[:12][::-1], cp_validation_exp[:12][::-1])
                    cp_validation_exp_interpolated = np.concatenate((cp_interpolated_exp_inf, cp_interpolated_exp))
                    sub_ax.scatter(x_airfoil_normalized_full, cp_validation_exp_interpolated, marker="^", label="Validation experimental", s=5)
                    # original
                    # sub_ax.scatter(x_validation_exp, cp_validation_exp, marker="^", label="Validation experimental orig", s=5)

                #### POTENTIAL ######
                potential_skin_data_filename = f'../reference_data/onera/potential_solver/cp_{labels[idx]}_new.dat'
                if os.path.exists(potential_skin_data_filename) and POTENTIAL:
                    x_validation_pot  = np.loadtxt(potential_skin_data_filename, usecols=(0,))
                    cp_validation_pot = np.loadtxt(potential_skin_data_filename, usecols=(1,))
                    cp_interpolated_pot = np.interp(x_airfoil_normalized, x_validation_pot[potential_lines[idx]+1:], cp_validation_pot[potential_lines[idx]+1:])
                    cp_interpolated_pot_inf = np.interp(x_airfoil_normalized, x_validation_pot[:potential_lines[idx]], cp_validation_pot[:potential_lines[idx]])
                    cp_validation_pot_interpolated = np.concatenate((cp_interpolated_pot_inf, cp_interpolated_pot))
                    sub_ax.scatter(x_airfoil_normalized_full, cp_validation_pot_interpolated, marker="*", label="Validation potential solver", s=5)
                    # original
                    # sub_ax.scatter(x_validation_pot, cp_validation_pot, marker="*", label="Validation potential solver orig", s=5)

            #### FOM ######
            if os.path.exists(fom_skin_data_filename) and FOM:
                # Interpola los valores de x y cp para la posición y_target
                cp_interpolated_fom = griddata((x_fom, y_fom), cp_fom, (x_grid, y_target), method='linear', fill_value=0.0)
                cp_interpolated_fom_inf = griddata((x_fom_inf, y_fom_inf), cp_fom_inf, (x_grid, y_target), method='linear')
                cp_interpolated_fom_full = np.concatenate((cp_interpolated_fom_inf, cp_interpolated_fom))
                if capture_directory == 'Validation':
                    error_exp_fom = np.linalg.norm(cp_interpolated_fom_full-cp_validation_exp_interpolated)/np.linalg.norm(cp_validation_exp_interpolated)
                    sub_ax.scatter(x_airfoil_normalized_full, -cp_interpolated_fom_full, marker="s", label=f"FOM-Exp e: {error_exp_fom:.2E}", s=5)

                    error_pot_fom = np.linalg.norm(cp_interpolated_fom_full-cp_validation_pot_interpolated)/np.linalg.norm(cp_validation_pot_interpolated)
                    sub_ax.scatter(x_airfoil_normalized_full, -cp_interpolated_fom_full, marker="s", label=f"FOM-Pot e: {error_pot_fom:.2E}", s=5)
                else:
                    sub_ax.scatter(x_airfoil_normalized_full, -cp_interpolated_fom_full, marker="s", label=f"FOM", s=5)

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

            sub_ax.set_title(f'Section y = {y_target}')
            # Ajustar las escalas de los ejes
            sub_ax.set_xlabel('x')
            sub_ax.set_ylabel('-Cp')
            if capture_directory == 'Validation':
                plt.legend(bbox_to_anchor=(-0.04, 4.8, 1.0, 0.1), loc='upper left', borderaxespad=0.0)
            else:
                plt.legend(bbox_to_anchor=(0.01, 4.8, 1.0, 0.1), loc='upper left', borderaxespad=0.0)

        # Mostrar el gráfico
        # plt.show()
        plt.savefig(capture_filename, dpi=400)
        plt.close('all')


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


if __name__ == "__main__":

    mu_train, mu_test = load_mu_parameters()

    Plot_Cps(mu_train, 'Train_Captures')
    Plot_Cps(mu_test, 'Test_Captures')

    mu = []
    mu.append([3.06,0.839])

    Plot_Cps(mu, 'Validation')
