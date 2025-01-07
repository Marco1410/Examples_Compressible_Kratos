import os
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import numpy as np


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

mu_train = []
mu_test  = []

if os.path.exists("mu_train_errors.npy"):
    mu_train = np.load('mu_train_errors.npy')
    mu_train =  [mu.tolist() for mu in mu_train]
if os.path.exists("mu_test_errors.npy"):
    mu_test = np.load('mu_test_errors.npy')
    mu_test =  [mu.tolist() for mu in mu_test]

markercolor = ["ob","xr","+g","sc","*m","Dy"]
tolerancias = [1e-3, 1e-5]
fig = plt.figure()
fig.set_figwidth(12.0)
fig.set_figheight(8.0)
if len(mu_train)>0:
    items_filtrados_1 = [item for item in mu_train if item[2] >= tolerancias[0]]
    items_filtrados_11 = [item for item in mu_train if item[2] < tolerancias[0]]
    items_filtrados_2 = [item for item in items_filtrados_11 if item[2] >= tolerancias[1]]
    items_filtrados_3 = [item for item in mu_train if item[2] < tolerancias[1]]

    train_a1 = [item[0] for item in items_filtrados_1]
    train_a2 = [item[0] for item in items_filtrados_2]
    train_a3 = [item[0] for item in items_filtrados_3]

    train_m1 = [item[1] for item in items_filtrados_1]
    train_m2 = [item[1] for item in items_filtrados_2]
    train_m3 = [item[1] for item in items_filtrados_3]

    if len(train_a1) > 0:
        fig = plt.plot(train_m1, train_a1, markercolor[0], label=f"train error > {tolerancias[0]}")
    if len(train_a2) > 0:
        fig = plt.plot(train_m2, train_a2, markercolor[1], label=f"train {tolerancias[0]} > error > {tolerancias[1]}")
    if len(train_a3) > 0:
        fig = plt.plot(train_m3, train_a3, markercolor[2], label=f"train error < {tolerancias[1]}")
    
    x = np.linspace(0.7, 0.75, 25)
    fig = plt.plot(x, -28.889*x+21.522, label=f"Limit")

if len(mu_test)>0:
    items_filtrados_1 = [item for item in mu_test if item[2] >= tolerancias[0]]
    items_filtrados_11 = [item for item in mu_test if item[2] < tolerancias[0]]
    items_filtrados_2 = [item for item in items_filtrados_11 if item[2] >= tolerancias[1]]
    items_filtrados_3 = [item for item in mu_test if item[2] < tolerancias[1]]

    test_a1 = [item[0] for item in items_filtrados_1]
    test_a2 = [item[0] for item in items_filtrados_2]
    test_a3 = [item[0] for item in items_filtrados_3]

    test_m1 = [item[1] for item in items_filtrados_1]
    test_m2 = [item[1] for item in items_filtrados_2]
    test_m3 = [item[1] for item in items_filtrados_3]

    if len(test_a1) > 0:
        fig = plt.plot(test_m1, test_a1, markercolor[3], label=f"test error > {tolerancias[0]}")
    if len(test_a2) > 0:
        fig = plt.plot(test_m2, test_a2, markercolor[4], label=f"test {tolerancias[0]} > error > {tolerancias[1]}")
    if len(test_a3) > 0:
        fig = plt.plot(test_m3, test_a3, markercolor[5], label=f"test error < {tolerancias[1]}")

fig = plt.title('Mu Values and approximation errors FOM vs ROM')
fig = plt.ylabel('Alpha')
fig = plt.xlabel('Mach')
fig = plt.grid(True)
fig = plt.legend(bbox_to_anchor=(.85, 1.03, 1., .102), loc='upper left', borderaxespad=0.)
fig = plt.savefig("MuValuesWithErrors.png")
fig = plt.close('all')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #