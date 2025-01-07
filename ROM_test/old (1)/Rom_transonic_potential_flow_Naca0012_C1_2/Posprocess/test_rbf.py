import os
import numpy as np
import matplotlib.pyplot as plt

class Rbf:
    def __init__(self, *args, epsilon=1.0):
        """
        Inicializa el interpolador RBF.
        args: Lista de arrays que contienen las coordenadas de cada punto (ángulo, Mach, etc.).
        epsilon: Parámetro de forma para la función gaussiana.
        """
        # Puntos conocidos
        self.x = np.vstack(args[:-1]).T  # Coordenadas de entrada (ángulo, mach, etc.)
        self.y = args[-1]  # Valores de los vectores conocidos
        self.epsilon = epsilon
        
        # Construir la matriz A
        self.A = self._compute_matrix(self.x)
        
        # Resolver para los pesos w
        self.w = np.linalg.solve(self.A, self.y)
    
    def _compute_matrix(self, x):
        """
        Calcula la matriz A utilizando la función gaussiana.
        """
        A = np.zeros((x.shape[0], x.shape[0]))
        for i in range(x.shape[0]):
            for j in range(x.shape[0]):
                A[i, j] = self._gaussian_kernel(np.linalg.norm(x[i] - x[j]))
        return A
    
    def _gaussian_kernel(self, r):
        """
        Función base radial Gaussiana.
        """
        return np.exp(-(self.epsilon * r) ** 2)
    
    def __call__(self, *args):
        """
        Interpola el valor en un nuevo punto.
        args: Las coordenadas del nuevo punto (nuevo ángulo, nuevo Mach, etc.).
        """
        x_new = np.vstack(args).T
        
        interpolated_values = []
        
        for x_pt in x_new:
            # Calcular el valor interpolado como combinación lineal de los valores base
            interpolated_value = 0
            for i in range(self.x.shape[0]):
                r = np.linalg.norm(self.x[i] - x_pt)
                interpolated_value += self.w[i] * self._gaussian_kernel(r)
            interpolated_values.append(interpolated_value)
        
        return np.array(interpolated_values)

def cargar_vectores_y_parametros(directorio):
    vectores = []
    parametros = []
    
    # Iterar sobre los archivos en el directorio
    for archivo in os.listdir(directorio):
        if archivo.endswith('.dat'):
            # Extraer ángulo y Mach del nombre del archivo
            nombre, _ = os.path.splitext(archivo)
            angulo, mach = map(float, nombre.split(', '))
            parametros.append([angulo, mach])
            
            # Cargar el vector del archivo .npy
            vector = np.loadtxt(os.path.join(directorio, archivo), usecols=(3,))
            vectores.append(vector)
    
    return np.array(parametros), np.array(vectores), np.loadtxt(os.path.join(directorio, archivo), usecols=(0,))

def interpolar_rbf_multidimensional(parametros, vectores, nuevo_angulo, nuevo_mach, epsilon=1.0):
    # Extraer ángulos y Machs
    angulos = parametros[:, 0]
    machs = parametros[:, 1]
    
    # Interpolar cada componente del vector por separado
    vectores_interpolados = []
    
    for i in range(vectores.shape[1]):  # Asumiendo que vectores tiene la forma (num_datos, dimension_vector)
        # Crear un interpolador RBF para cada componente del vector
        rbf = Rbf(angulos, machs, vectores[:, i], epsilon=epsilon)
        # Interpolar para los nuevos valores de ángulo y Mach
        interpolado_componente = rbf(nuevo_angulo, nuevo_mach)
        vectores_interpolados.append(interpolado_componente)
    
    return np.array(vectores_interpolados)

# Directorio donde se almacenan los archivos .npy
directorio = 'FOM_Skin_Data'

# Cargar parámetros (ángulo, mach) y los vectores correspondientes
parametros, vectores, x = cargar_vectores_y_parametros(directorio)

# Definir nuevos valores de ángulo y Mach para la interpolación
nuevo_angulo = 1.2
nuevo_mach = 0.73

# Obtener el vector interpolado
vector_interpolado = interpolar_rbf_multidimensional(parametros, vectores, nuevo_angulo, nuevo_mach, epsilon=0.2)

# Cp = np.load(f'rom_data/RightBasisMatrix.npy').dot(vector_interpolado)

Cp = vector_interpolado

plt.plot(x,-Cp,'bs',label="Cp")
plt.title('RBF Test')
plt.ylabel('Cp')
plt.xlabel('x')
plt.grid(True)
plt.legend()
plt.show()
plt.close('all')