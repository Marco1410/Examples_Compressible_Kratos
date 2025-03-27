from parameters_manager import *  
from database_rbf import *



if __name__ == "__main__":

    # Posprocess -- Plot Cp distributions 2d
    
    weight          = 1.0     
    alpha           = 0.85   
    mu_id           = 0   
    svd_tol         = 1e-12  
    solver_strategy = 1     
    mu              = [1.0, 0.75]
    case            = 'FOM_Run'

    target_directory = 'output_plots'
    os.makedirs(target_directory, exist_ok=True)

    case_name = f'{target_directory}/2D_{case}_{mu[0]}, {mu[1]}.png'
    
    vtk_filename = f"Results/vtk_output_{case}{mu[0]}, {mu[1]}, {weight}, {alpha}, {mu_id}, {svd_tol}, {solver_strategy}/MainModelPart_0_1.vtk" 

    mesh = pv.read(vtk_filename).extract_geometry()
    points, cells = mesh.points, mesh.faces.reshape(-1, 4)[:, 1:4]
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
    ax.set_title(f'Alpha={mu[0]:.3f}, Mach={mu[1]:.3f}', size=13, pad=10)
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(-0.75, 1.0)
    plt.savefig(case_name, dpi=200)
    plt.close('all')


    # # Posprocess -- Plot ACTIVATION_LEVEL distribution 2D (Elemental)
    
    weight          = 1.0     
    alpha           = 0.85   
    mu_id           = 0   
    svd_tol         = 1e-12  
    solver_strategy = 1     
    mu              = [1.0, 0.75]
    case            = 'FOM_Run'

    target_directory = 'output_plots'
    os.makedirs(target_directory, exist_ok=True)

    case_name = f'{target_directory}/2D_{case}_{mu[0]}, {mu[1]}_ACTIVATION_LEVEL.png'
    
    vtk_filename = f"Results/vtk_output_{case}{mu[0]}, {mu[1]}, {weight}, {alpha}, {mu_id}, {svd_tol}, {solver_strategy}/MainModelPart_0_1.vtk" 

    mesh = pv.read(vtk_filename).extract_geometry()
    points, cells = mesh.points, mesh.faces.reshape(-1, 4)[:, 1:4]

    # Extraer la variable ACTIVATION_LEVEL desde los datos elementales
    activation_level = mesh.cell_data['ACTIVATION_LEVEL']

    # Asignar valores elementales a los nodos promediando valores de los elementos vecinos
    activation_nodal = np.zeros(len(points))
    count = np.zeros(len(points))

    for i, cell in enumerate(cells):
        activation_nodal[cell] += activation_level[i]
        count[cell] += 1

    activation_nodal /= np.maximum(count, 1)  # Evitar división por cero

    # fig, ax = plt.subplots(figsize=(10, 8))
    # triang = tri.Triangulation(points[:, 0], points[:, 1], cells)

    # cmap = cm.get_cmap('coolwarm')  # viridis coolwarm
    
    # # contour = ax.tricontour(triang, activation_nodal, levels=15, cmap='coolwarm', linewidths=0.8)
    # ax.tricontourf(triang, activation_nodal, levels=100, cmap=cmap, alpha=1.0)
    
    # # cbar = plt.colorbar(cm.ScalarMappable(norm=plt.Normalize(vmin=activation_nodal.min(), vmax=activation_nodal.max()), cmap=cmap), ax=ax, label='ACTIVATION LEVEL')
    # # cbar.ax.tick_params(labelsize=10)
    
    # # ax.set_xlabel('X')
    # # ax.set_ylabel('Y')
    # # ax.set_title(f'Alpha={mu[0]:.3f}, Mach={mu[1]:.3f}', size=13, pad=10)
    # ax.set_xlim(-0.5, 1.5)
    # ax.set_ylim(-0.75, 1.0)
    
    # plt.savefig(case_name, dpi=200)
    # plt.close('all')

    fig, ax = plt.subplots(figsize=(10, 8))
    triang = tri.Triangulation(points[:, 0], points[:, 1], cells)

    cmap = cm.get_cmap('coolwarm')  # Color map

    # Dibujar la distribución de ACTIVATION_LEVEL
    ax.tricontourf(triang, activation_nodal, levels=100, cmap=cmap, alpha=1.0)

    # **Dibujar la malla en negro**
    ax.triplot(triang, 'k-', lw=0.5)  # 'k-' = color negro, lw=0.5 = línea delgada

    # Configuración del gráfico
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(-0.75, 1.0)

    plt.savefig(case_name, dpi=200)
    plt.close('all')

