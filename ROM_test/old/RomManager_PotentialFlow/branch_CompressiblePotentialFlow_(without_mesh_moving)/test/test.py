import numpy as np

# Read the CSV file
phi_matrix = np.genfromtxt('eigen_mPhiGlobal.csv', delimiter=',')

print(" ")
print(" phi ")
print("shape:",phi_matrix.shape)
print("norm:",np.linalg.norm(phi_matrix))

u_matrix = np.load("u_matrix.npy")

filas_no_cero = np.any(u_matrix != 0, axis=1)
numero_filas_no_cero = np.sum(filas_no_cero)

print(" ")
print(" u ")
print("shape:",u_matrix.shape)
print("norm:",np.linalg.norm(u_matrix))
print("non cero rows:", numero_filas_no_cero)

u_matrix_filtered = u_matrix[u_matrix[:, 0] != 0]

print(" ")
print(" u filtered ")
print("shape:",u_matrix_filtered.shape)
print("norm:",np.linalg.norm(u_matrix_filtered))

S_matrix = np.load("s_matrix.npy")

print(" ")
print(" s ")
print("shape:",S_matrix.shape)

V_matrix = np.load("v_matrix.npy")

print(" ")
print(" v ")
print("shape:",V_matrix.shape)

snapshotsmatrix = np.load("snapshotsmatrix.npy")

print(" ")
print(" sanpshotsmatrix ")
print("shape:",snapshotsmatrix.shape)

print(" ")
print(" norm(Snapshotsmatrix - U*S*Vt) = ", np.linalg.norm(snapshotsmatrix - np.dot(u_matrix * S_matrix , V_matrix.transpose())))

aux = np.dot(u_matrix.transpose(),snapshotsmatrix)
print(" (norm(Sn - U*Ut*Sn)) / norm(Sn) = ", (np.linalg.norm(snapshotsmatrix - np.dot(u_matrix,aux)))/np.linalg.norm(snapshotsmatrix))

aux = np.dot(phi_matrix.transpose(),snapshotsmatrix)
print(" (norm(Sn - phi*phit*Sn)) / norm(Sn) = ", (np.linalg.norm(snapshotsmatrix - np.dot(phi_matrix,aux)))/np.linalg.norm(snapshotsmatrix))