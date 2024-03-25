import numpy as np

# Definir los centroides (puedes añadir más o cambiar estos valores)
centroides = np.array([(10, 30), (20, 30), (40, 50)])

# Definir el punto de interés
punto = (25, 35)

# Calcular la distancia euclidiana desde el punto a cada centroide
distancias = np.sqrt(np.sum((centroides - punto) ** 2, axis=1))

# Encontrar la distancia mínima y el centroide más cercano
min_distancia = np.min(distancias)
id_min_distancia = np.argmin(distancias)
centroide_mas_cercano = centroides[np.argmin(distancias)]

print(f"Distancias: {distancias}")
print(f"distancia 1: {np.linalg.norm(centroides[0]-punto)}")
print(f"distancia 2: {np.linalg.norm(centroides[1]-punto)}")
print(f"distancia 3: {np.linalg.norm(centroides[2]-punto)}")
print(f"Centroide más cercano: {centroide_mas_cercano}")
print(f"Distancia mínima: {min_distancia}")
print(f"Id dentro del vector: {id_min_distancia}")