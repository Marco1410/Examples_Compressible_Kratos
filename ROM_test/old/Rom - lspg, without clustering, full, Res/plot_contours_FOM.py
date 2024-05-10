import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os

# Función para cargar datos de un archivo específico.
def load_data(filename):
    # Cargamos las columnas correspondientes a 'x', 'y' y 'Cp' del archivo .dat.
    return np.loadtxt(filename, usecols=(0, 1, 3))

# Función para generar el gráfico
def plot_data(data_folder):
    # Lista para almacenar los datos para y > 0 e y < 0
    datasets_y_positive = []
    datasets_y_negative = []
    angles = []  # Lista para guardar los valores de los ángulos de ataque

    # Itera a través de los archivos en la carpeta de datos
    for file in os.listdir(data_folder):
        if file.endswith('.dat'):
            filepath = os.path.join(data_folder, file)
            data = load_data(filepath)
            # Ordenar los datos por la columna 'x'
            indices_sorted = np.argsort(data[:, 0])
            data_sorted = data[indices_sorted]

            # Separamos los datos en y > 0 y y < 0
            y_positive = data_sorted[data_sorted[:, 1] > 0]
            y_negative = data_sorted[data_sorted[:, 1] < 0]

            if y_positive.size > 0:
                datasets_y_positive.append(y_positive)
            if y_negative.size > 0:
                datasets_y_negative.append(y_negative)
            
            # Extraemos el ángulo de ataque del nombre del archivo y lo convertimos a float
            angle, _ = os.path.splitext(file)[0].split(',')
            angles.append(float(angle.strip()))

    # Obtener los colores para cada ángulo de ataque basado en una escala
    angles_sorted_indices = np.argsort(angles)
    colors = cm.viridis(np.linspace(0, 1, len(angles)))

    # Crear los subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

    # Plot para y > 0
    for i, data in enumerate(np.array(datasets_y_positive)[angles_sorted_indices]):
        x = data[:, 0]  # Posición x
        cp = data[:, 2]  # Coeficiente de presión Cp
        label = f"{angles[angles_sorted_indices[i]]:.2f}°"
        ax1.plot(x, cp, label=label, color=colors[i])

    # Plot para y < 0
    for i, data in enumerate(np.array(datasets_y_negative)[angles_sorted_indices]):
        x = data[:, 0]  # Posición x
        cp = data[:, 2]  # Coeficiente de presión Cp
        label = f"{angles[angles_sorted_indices[i]]:.2f}°"
        ax2.plot(x, cp, label=label, color=colors[i])

    # Configuraciones de los subplots para y > 0
    ax1.set_xlabel('x')
    ax1.set_ylabel('Cp')
    ax1.set_title('Cp vs x for y > 0')
    # ax1.legend(title="Ángulo de ataque")
    ax1.grid(True)
    ax1.invert_yaxis()

    # Configuraciones de los subplots para y < 0
    ax2.set_xlabel('x')
    ax2.set_ylabel('Cp')
    ax2.set_title('Cp vs x for y < 0')
    # ax2.legend(title="Ángulo de ataque")
    ax2.grid(True)
    ax2.invert_yaxis()

    # Ajuste del layout y mostrar la gráfica
    plt.tight_layout()
    plt.show()

# Llamar a la función de trazado con la carpeta que contiene los datos
data_folder = 'DataBase/Data'
plot_data(data_folder)
