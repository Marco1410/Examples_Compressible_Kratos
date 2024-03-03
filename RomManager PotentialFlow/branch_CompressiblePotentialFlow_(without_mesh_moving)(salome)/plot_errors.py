import os
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import pickle
import numpy as np


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# load mu parameters
#

def load_mu_parameters(with_errors=False):
    if with_errors:
        if os.path.exists("mu_train_errors.dat") and os.path.exists("mu_test_errors.dat"):
            archivo = open('mu_train_errors.dat', 'rb')
            mu_train = pickle.load(archivo)
            archivo.close()
            archivo = open('mu_test_errors.dat', 'rb')
            mu_test = pickle.load(archivo)
            archivo.close()
        elif os.path.exists("mu_train_errors.dat"):
            archivo = open('mu_train_errors.dat', 'rb')
            mu_train = pickle.load(archivo)
            archivo.close()
            mu_test = []
        elif os.path.exists("mu_test_errors.dat"):
            archivo = open('mu_test_errors.dat', 'rb')
            mu_test = pickle.load(archivo)
            archivo.close()
            mu_train = []
    else:    
        if os.path.exists("mu_train.dat") and os.path.exists("mu_test.dat"):
            archivo = open('mu_train.dat', 'rb')
            mu_train = pickle.load(archivo)
            archivo.close()
            archivo = open('mu_test.dat', 'rb')
            mu_test = pickle.load(archivo)
            archivo.close()
        elif os.path.exists("mu_train.dat"):
            archivo = open('mu_train.dat', 'rb')
            mu_train = pickle.load(archivo)
            archivo.close()
            mu_test = []
        elif os.path.exists("mu_test.dat"):
            archivo = open('mu_test.dat', 'rb')
            mu_test = pickle.load(archivo)
            archivo.close()
            mu_train = []
    return mu_train, mu_test

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# plot Cps
#

print("Loading parameters")
mu_train, mu_test = load_mu_parameters(with_errors=True)
print(" ")
print(len(mu_train), " mu train parameter loaded")
if len(mu_train) > 0: print(mu_train)
print(" ")
print(len(mu_test), " mu test parameter loaded")
if len(mu_test) > 0: print(mu_test)
print(" ")

markercolor = ["ob","xr","+g","sc","*m","Dy"]
tolerancias = [1.0, 1e-3]
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
        fig = plt.plot(train_m1, train_a1, markercolor[0], label="train error > 1%")
    if len(train_a2) > 0:
        fig = plt.plot(train_m2, train_a2, markercolor[1], label="train 1% > error > 1e-3%")
    if len(train_a3) > 0:
        fig = plt.plot(train_m3, train_a3, markercolor[2], label="train error < 1e-3%")

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
        fig = plt.plot(test_m1, test_a1, markercolor[3], label="test error > 1%")
    if len(test_a2) > 0:
        fig = plt.plot(test_m2, test_a2, markercolor[4], label="test 1% > error > 1e-3%")
    if len(test_a3) > 0:
        fig = plt.plot(test_m3, test_a3, markercolor[5], label="test error < 1e-3%")

fig = plt.title('Mu Values and approximation errors FOM vs ROM')
fig = plt.ylabel('Alpha')
fig = plt.xlabel('Mach')
fig = plt.grid(True)
fig = plt.legend(bbox_to_anchor=(.85, 1.03, 1., .102), loc='upper left', borderaxespad=0.)
fig = plt.savefig("MuValuesWithErrors.png")
fig = plt.close('all')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #