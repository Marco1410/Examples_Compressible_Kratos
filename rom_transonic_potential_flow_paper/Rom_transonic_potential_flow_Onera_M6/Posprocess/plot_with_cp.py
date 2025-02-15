import os
import numpy as np
import pyvista as pv 
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# load parameters
#
def load_mu_parameters(name):
    filename = f'{name}.npy'
    if os.path.exists(filename):
        mu_npy = np.load(filename)
        mu =  [mu.tolist() for mu in mu_npy]
    else:
        mu = []
    return mu

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



if __name__ == "__main__":

    strategies = ['galerkin','lspg']
    n = [0, 0]
    indices_list = [[1,9,4,6], [1,9,8,2]]

    for id, strategy, indices in zip(n, strategies, indices_list):
        mu_train      = load_mu_parameters(f'Mu_history/{id}_{strategy}_mu_train')
        mu_test       = load_mu_parameters(f'Mu_history/{id}_{strategy}_mu_test')
        mu_validation = load_mu_parameters(f'Mu_history/mu_validation')

        mach_train = [mu[1] for mu in mu_train]
        alpha_train = [mu[0] for mu in mu_train]
        mach_test = [mu[1] for mu in mu_test]
        alpha_test = [mu[0] for mu in mu_test]
        mach_validation = [mu[1] for mu in mu_validation]
        alpha_validation = [mu[0] for mu in mu_validation]

        # Dibujar subplots
        # Crear figura y ejes
        fig = plt.figure(figsize=(10, 10))

        # Posición y tamaño del gráfico central (ajustable)
        central_plot_position = [0.35, 0.35, 0.3, 0.3]  # [left, bottom, width, height]
        ax = fig.add_axes(central_plot_position)
        ax.scatter(mach_train, alpha_train, color='blue', label='Train')
        ax.scatter(mach_test, alpha_test, color='red', label='Test')
        ax.scatter(mach_validation, alpha_validation, color='green', label='Validation')
        ax.set_xlim(0.8,0.85)
        ax.set_xlabel('Mach')
        ax.set_ylabel('Alpha')
        ax.legend(bbox_to_anchor=(1.01, 0.9, 1.0, 0.1), loc='upper left', borderaxespad=0.0)
        ax.set_title('Parameter space', size=13, pad=8)

        # Posiciones relativas de los subplots
        positions = [
            [0.05, 0.70, 0.25, 0.25],
            [0.68, 0.70, 0.25, 0.25],
            [0.05, 0.05, 0.25, 0.25],
            [0.68, 0.05, 0.25, 0.25]
        ]

        # Dibujar subplots
        I=II=III=IV=False
        for idx in indices:
            filename = f'FOM_Skin_Data/{alpha_train[idx]}, {mach_train[idx]}.dat'
            if (mach_train[idx] <= np.min(mach_train) + ((np.max(mach_train)-np.min(mach_train))/2)):
                if (alpha_train[idx] <= np.min(alpha_train) + ((np.max(alpha_train)-np.min(alpha_train))/2)):
                    if not III:
                        pos = positions[2]
                        III = True
                    else: 
                        pos = positions[0]
                else:
                    if not I:
                        pos = positions[0]
                        I = True
                    else: 
                        pos = positions[2]
            else:
                if (alpha_train[idx] <= np.min(alpha_train) + ((np.max(alpha_train)-np.min(alpha_train))/2)):
                    if not IV:
                        pos = positions[3]
                        IV = True
                    else: 
                        pos = positions[1]
                else:
                    if not II:
                        pos = positions[1]
                        II = True
                    else: 
                        pos = positions[3]

            sub_ax = fig.add_axes(pos, projection='3d')
            x_fom = np.loadtxt(filename, usecols=(0,))
            y_fom = np.loadtxt(filename, usecols=(1,))
            z_fom = np.loadtxt(filename, usecols=(2,))
            cp_fom = np.loadtxt(filename, usecols=(3,))

            # Cargar el archivo .vtk
            mesh = pv.read(f"Results/vtk_output_FOM_Fit{idx}/MainModelPart_Wing_0_1.vtk")
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
            sub_ax.add_collection3d(tri_collection)
            sub_ax.view_init(elev=30, azim=135) 

            sub_ax.set_title(f'Angle={np.round(alpha_train[idx],3)}, Mach={np.round(mach_train[idx],3)}')
            # Ajustar las escalas de los ejes
            sub_ax.set_zlim([-0.5, 0.5])
            sub_ax.set_ylim([ 0.0, 1.5])
            sub_ax.set_xlabel('x')
            sub_ax.set_ylabel('y')
            sub_ax.set_zlabel('z')
            
            # Conectar con líneas o flechas
            arrowprops = dict(arrowstyle="->", color='grey', lw=1)
            ax.annotate('',
                        xy=(mach_train[idx], alpha_train[idx]), 
                        xycoords='data',
                        xytext=(pos[0] + pos[2]/2, pos[1] + pos[3]/2), 
                        textcoords='figure fraction',
                        arrowprops=arrowprops)

        # Mostrar el gráfico
        # plt.show()
        plt.savefig(f'Mu_parameters_cp_{strategy}.pdf', dpi=400)
