import meshio

mesh_mdpa = meshio.read("naca0012_0aoa_00.mdpa")

# Resumen de la malla MDPA
print(f"Puntos MDPA: {mesh_mdpa.points.shape}")
print(f"Tipos de células MDPA: {set(cells.type for cells in mesh_mdpa.cells)}")

mesh_mdpa.write("test","med")


mesh_med = meshio.read("test.med")

# Resumen de la malla med
print(f"Puntos med: {mesh_med.points.shape}")
print(f"Tipos de células med: {set(cells.type for cells in mesh_med.cells)}")

mesh_med.write("test","mdpa")