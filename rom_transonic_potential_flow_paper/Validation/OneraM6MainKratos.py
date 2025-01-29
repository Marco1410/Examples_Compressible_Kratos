#####################################################################################################
#####################################################################################################
######################                ONERA M6 VALIDATION
#####################################################################################################
###################### This data is from:
######################  
###################### Experimental:
######################     V. Schmitt and F. Charpin. Pressure Distributions on the
######################     ONERA-M6-Wing at Transonic Mach Numbers, Experimental Data
######################     Base for Computer Program Assessment. Tech. rep. Report of the
######################     Fluid Dynamics Panel Working Group 04, AGARD AR 138, 1979
######################     https://www.sto.nato.int/publications/AGARD/AGARD-AR-138/AGARD-AR-138.pdf
######################
###################### FV Potential solution:
######################     @article{LYU2017951,
######################     title = {A fast and automatic full-potential finite volume solver 
######################              on Cartesian grids for unconventional configurations},
######################     doi = {https://doi.org/10.1016/j.cja.2017.03.001},
######################     url = {https://www.sciencedirect.com/science/article/pii/S1000936117300730},
######################     author = {Fanxi LYU and Tianhang XIAO and Xiongqing YU},
######################     }
######################
###################### FE Potential solution:
######################     @phdthesis{dissertation,
######################         author = {López Canalejo, Iñigo Pablo},
######################         title = {A finite-element transonic potential flow solver with an 
######################                  embedded wake approach for aircraft conceptual design},
######################         url = {https://mediatum.ub.tum.de/1633175},
######################     }
#####################################################################################################

import os
import sys
import time
import importlib
import numpy as np
import pyvista as pv 
import matplotlib.cm as cm
import matplotlib.tri as tri
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import KratosMultiphysics
import KratosMultiphysics.kratos_utilities
import KratosMultiphysics.RomApplication
import KratosMultiphysics.FluidDynamicsApplication 


def CreateAnalysisStageWithFlushInstance(cls, global_model, parameters):
    class AnalysisStageWithFlush(cls):

        def __init__(self, model,project_parameters, flush_frequency=10.0):
            super().__init__(model,project_parameters)
            self.flush_frequency = flush_frequency
            self.last_flush = time.time()
            sys.stdout.flush()

        def Initialize(self):
            super().Initialize()
            sys.stdout.flush()

        def FinalizeSolutionStep(self):
            super().FinalizeSolutionStep()
    
            if self.parallel_type == "OpenMP":
                now = time.time()
                if now - self.last_flush > self.flush_frequency:
                    sys.stdout.flush()
                    self.last_flush = now

    return AnalysisStageWithFlush(global_model, parameters)


if __name__ == "__main__":

    UPDATE_DATA     = True
    PLOTCPS         = True
    PLOT2D          = True
    PLOT3D          = True
    ##########################
    # SELECT WHICH MODEL TO SHOW
    FE_KRATOS       = True
    FE_KRATOS_INIGO = False
    EXPERIMENTAL    = True
    FV_POTENTIAL    = True
    ##########################

    mu_validation = []
    mu_validation.append([3.06,0.839])

    with open("OneraM6ProjectParameters.json", 'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    analysis_stage_module_name = parameters["analysis_stage"].GetString()
    analysis_stage_class_name = analysis_stage_module_name.split('.')[-1]
    analysis_stage_class_name = ''.join(x.title() for x in analysis_stage_class_name.split('_'))

    analysis_stage_module = importlib.import_module(analysis_stage_module_name)
    analysis_stage_class = getattr(analysis_stage_module, analysis_stage_class_name)

    for mu in mu_validation:

        angle_of_attack = mu[0]
        mach_infinity   = mu[1]

        case_name = f'{angle_of_attack}, {mach_infinity}'

        wake_normal     = [-np.sin(angle_of_attack*np.pi/180),0.0,np.cos(angle_of_attack*np.pi/180)]
        parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(np.double(angle_of_attack))
        parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(np.double(mach_infinity))
        parameters["processes"]["boundary_conditions_process_list"][1]["Parameters"]["wake_process_cpp_parameters"]["wake_normal"].SetVector(wake_normal)
        parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString(f'Results_OneraM6/{angle_of_attack}, {mach_infinity}')
        parameters["output_processes"]["vtk_output"][0]["Parameters"]["output_path"].SetString(f'Results_OneraM6/{angle_of_attack}, {mach_infinity}')

        kratos_skin_data_filename = f"Onera_M6_kratos_Cp_data.dat"
        if not os.path.exists(kratos_skin_data_filename) or UPDATE_DATA:
            global_model = KratosMultiphysics.Model()
            simulation = CreateAnalysisStageWithFlushInstance(analysis_stage_class, global_model, parameters)
            simulation.Run()
        else:
            x  = np.loadtxt(kratos_skin_data_filename, usecols=(0,))
            y  = np.loadtxt(kratos_skin_data_filename, usecols=(1,))
            z  = np.loadtxt(kratos_skin_data_filename, usecols=(2,))
            cp = np.loadtxt(kratos_skin_data_filename, usecols=(3,))

        # LOAD .vtk FILE
        mesh = pv.read(f"Results_OneraM6/{angle_of_attack}, {mach_infinity}/MainModelPart_Wing_0_1.vtk")

        # EXTRACT SURFACE GEOMETRY
        surface_mesh = mesh.extract_geometry()

        # EXTRACT POINTS AND CELLS
        points = surface_mesh.points
        cells = surface_mesh.faces.reshape(-1, 4)[:, 1:4]

        # EXTRACT 'PRESSURE_COEFFICIENT'
        pressure_coeff = surface_mesh.point_data['PRESSURE_COEFFICIENT']

        # CREATE COLOR MAP WITH PRESSURE COEFFICIENT
        norm = plt.Normalize(vmin=pressure_coeff.min(), vmax=pressure_coeff.max())
        cmap = cm.get_cmap('viridis')

        ########## 3D PLOT #####################################################################
        if PLOT3D:
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
            ax.set_title(f' Wing Pressure distribution - Angle={np.round(mu[0],3)}, Mach={np.round(mu[1],3)}', size=13, pad=8, y=-0.15)
            # plt.show()
            plt.savefig("Onera M6 Validation 3D view.png", dpi=200)
            plt.close('all')

        ########## 2D PLOT #####################################################################
        if PLOT2D:
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
            ax.set_title(f' Wing Pressure distribution - Angle={np.round(mu[0],3)}, Mach={np.round(mu[1],3)}', size=13, pad=8, y=-0.125)
            # plt.show()
            plt.savefig("Onera M6 Validation 2D view.png", dpi=200)
            plt.close('all')

        # WRITE CP SUBPLOTS 
        #########################################################
        if PLOTCPS:
            fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(13, 7))
            labels = [95, 90, 80, 65, 44, 20] 
            fv_potential_lines_in_data_files = [117, 116, 116, 116, 119, 122]
            fe_kratos_inigo_lines_in_data = [80, 73, 72, 74, 74, 81]

            if not os.path.exists(kratos_skin_data_filename) or UPDATE_DATA:
                modelpart = global_model["MainModelPart.Wing"]
                x  = np.zeros(modelpart.NumberOfNodes())
                y  = np.zeros(modelpart.NumberOfNodes())
                z  = np.zeros(modelpart.NumberOfNodes())
                cp = np.zeros(modelpart.NumberOfNodes())

                fout = open(kratos_skin_data_filename,'w')
                for i,node in enumerate(modelpart.Nodes):
                    x[i] = node.X; y[i] = node.Y; z[i] = node.Z
                    cp[i] = node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
                    fout.write("%s %s %s %s\n" %(x[i],y[i],z[i],cp[i]))
                fout.close()
                
            indexes    = z > 0
            x_k      = x[indexes]
            y_k      = y[indexes]
            cp_k     = cp[indexes]
            indexes    = z < 0
            x_k_inf  = x[indexes]
            y_k_inf  = y[indexes]
            cp_k_inf = cp[indexes]

            for idx, sub_ax in enumerate(axs.flat):

                # SET Y POSITION TO PLOT CP DISTRIBUTION
                b = 1.1963
                y_target = labels[idx]*b/100 # y_target = y/b

                # CREATE X LINEAR SPACE
                x_min = (y_target)*np.tan(30*np.pi/180)
                x_max = (y_target)*np.tan(15.8*np.pi/180)+0.8059
                x_grid = np.linspace(x_min, x_max, 80)
                # SET X BETWEEN 0 and 1
                x_airfoil_normalized = (x_grid - x_min) / (x_max - x_min)
                x_airfoil_normalized_full = np.concatenate((x_airfoil_normalized, x_airfoil_normalized))
                
                #### VALIDATION ######
                experimental_skin_data_filename = f'../reference_data/onera/experiment/cp_{labels[idx]}.dat'
                if os.path.exists(experimental_skin_data_filename) and EXPERIMENTAL:
                    x_validation_exp  = np.loadtxt(experimental_skin_data_filename, usecols=(0,))
                    cp_validation_exp = np.loadtxt(experimental_skin_data_filename, usecols=(1,))
                    # ORIGINAL
                    sub_ax.scatter(x_validation_exp, cp_validation_exp, marker="+", label="Experimental data", s=15, color='red')
                    ####################################################################################################################
                    # INTERPOLATED
                    cp_interpolated_exp = np.interp(x_airfoil_normalized, x_validation_exp[13:], cp_validation_exp[13:])
                    cp_interpolated_exp_inf = np.interp(x_airfoil_normalized, x_validation_exp[:12][::-1], cp_validation_exp[:12][::-1])
                    cp_validation_exp_interpolated = np.concatenate((cp_interpolated_exp_inf, cp_interpolated_exp))
                    # sub_ax.scatter(x_airfoil_normalized_full, cp_validation_exp_interpolated, marker="+", label="Experimental", s=7)

                #### FV POTENTIAL SOLVER ######
                potential_skin_data_filename = f'../reference_data/onera/potential_solver/cp_{labels[idx]}_new.dat'
                if os.path.exists(potential_skin_data_filename) and FV_POTENTIAL:
                    x_validation_pot  = np.loadtxt(potential_skin_data_filename, usecols=(0,))
                    cp_validation_pot = np.loadtxt(potential_skin_data_filename, usecols=(1,))
                    # ORIGINAL
                    sub_ax.scatter(x_validation_pot, cp_validation_pot, marker="*", label="FV Potential Solver", s=10)
                    ####################################################################################################################
                    # INTERPOLATED
                    cp_interpolated_pot = np.interp(x_airfoil_normalized, x_validation_pot[fv_potential_lines_in_data_files[idx]+1:], cp_validation_pot[fv_potential_lines_in_data_files[idx]+1:])
                    cp_interpolated_pot_inf = np.interp(x_airfoil_normalized, x_validation_pot[:fv_potential_lines_in_data_files[idx]], cp_validation_pot[:fv_potential_lines_in_data_files[idx]])
                    cp_validation_pot_interpolated = np.concatenate((cp_interpolated_pot_inf, cp_interpolated_pot))
                    # sub_ax.scatter(x_airfoil_normalized_full, cp_validation_pot_interpolated, marker="*", label="FV Potential solver", s=4)

                #### FE KRATOS IÑIGO ######
                potential_skin_data_filename = f'../reference_data/onera/kratos_potential_solver/kratos_cp_{labels[idx]}.dat'
                if os.path.exists(potential_skin_data_filename) and FE_KRATOS_INIGO:
                    x_validation_kratos  = np.loadtxt(potential_skin_data_filename, usecols=(0,))
                    cp_validation_kratos = np.loadtxt(potential_skin_data_filename, usecols=(1,))
                    # ORIGINAL
                    # sub_ax.scatter(x_validation_kratos, cp_validation_kratos, marker=".", label="FE kratos potential solver", s=4)
                    ####################################################################################################################
                    # INTERPOLATED
                    cp_interpolated_kratos = np.interp(x_airfoil_normalized, x_validation_kratos[fe_kratos_inigo_lines_in_data[idx]+1:], cp_validation_kratos[fe_kratos_inigo_lines_in_data[idx]+1:])
                    cp_interpolated_kratos_inf = np.interp(x_airfoil_normalized, x_validation_kratos[:fe_kratos_inigo_lines_in_data[idx]], cp_validation_kratos[:fe_kratos_inigo_lines_in_data[idx]])
                    cp_validation_kratos_interpolated = np.concatenate((cp_interpolated_kratos_inf, cp_interpolated_kratos))
                    sub_ax.scatter(x_airfoil_normalized_full, cp_validation_kratos_interpolated, marker=".", label="FE kratos Iñigo", s=4)

                #### FE KRATOS ######
                if FE_KRATOS:
                    # INTERPOLATED
                    cp_interpolated_fe_kratos = griddata((x_k, y_k), cp_k, (x_grid, y_target), method='linear', fill_value=0.25)
                    cp_interpolated_fe_kratos_inf = griddata((x_k_inf, y_k_inf), cp_k_inf, (x_grid, y_target), method='linear')
                    cp_interpolated_fe_kratos_full = np.concatenate((cp_interpolated_fe_kratos_inf, cp_interpolated_fe_kratos))
                    sub_ax.scatter(x_airfoil_normalized_full, -cp_interpolated_fe_kratos_full, marker="^", label=f"FE Kratos Potential Solver", s=7)

                # if FV_POTENTIAL and FE_KRATOS:
                #     error_pot_fe_kratos = np.linalg.norm(-cp_interpolated_fe_kratos_full-cp_validation_pot_interpolated)/np.linalg.norm(cp_validation_pot_interpolated)
                #     t = (f'Kratos vs FV pot. sol. \n Diff. {error_pot_fe_kratos:.2E}')
                #     sub_ax.text(0.25, -0.6, t, size=9, bbox=dict(boxstyle="round", facecolor='red', alpha=0.5))

                if idx == 2:
                    sub_ax.legend(bbox_to_anchor=(1.04, 1), loc='upper left', borderaxespad=0.0)
                sub_ax.set_title(f'Section y/b = {np.round(y_target/b,2)}')
                if idx == 3 or idx == 4 or idx == 5:
                    sub_ax.set_xlabel('x')
                if idx == 0 or idx == 3:
                    sub_ax.set_ylabel('-Cp')

            plt.tight_layout()
            # plt.show()
            plt.savefig("Onera M6 Validation 2D CPs view.png", dpi=200)
            plt.close('all')