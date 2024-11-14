import matplotlib.pyplot as plt

# Crear una figura y un eje
fig, ax = plt.subplots()

# Definir los rangos y los colores para cada región
regions = [
    ((0.70, 0.00), (0.73, 1.00), 'red'),    
    ((0.73, 0.00), (0.74, 1.00), 'green'),  
    ((0.74, 0.00), (0.75, 1.00), 'yellow'),  
    ((0.75, 0.00), (0.76, 1.00), 'magenta'),  
    ((0.71, 0.60), (0.73, 1.00), 'blue'),   
    ((0.70, 1.00), (0.73, 1.75), 'purple'),  
    ((0.73, 1.00), (0.74, 1.75), 'orange'),  
    ((0.74, 1.00), (0.75, 1.75), 'brown'),  
    ((0.75, 1.00), (0.76, 1.75), 'pink'),  
    ((0.70, 1.75), (0.73, 2.50), 'cyan'),  
    ((0.73, 1.75), (0.74, 2.50), 'gray'),  
    ((0.74, 1.75), (0.75, 2.50), 'lime'),  
    ((0.75, 1.75), (0.76, 2.50), 'olive')  
]

# Dibujar las regiones como rectángulos coloreados
for bottom_left, top_right, color in regions:
    rect = plt.Rectangle(bottom_left, top_right[0] - bottom_left[0], top_right[1] - bottom_left[1], 
                         facecolor=color, edgecolor='black', alpha=0.7)
    ax.add_patch(rect)

# Ajustar los límites de los ejes
ax.set_xlim(0.70, 0.76)
ax.set_ylim(0, 2.5)

# Etiquetas y título
plt.xlabel('Mach')
plt.ylabel('Alpha')
plt.title('Parameter setting regions')

# Mostrar la gráfica
plt.grid(True)
plt.show()

