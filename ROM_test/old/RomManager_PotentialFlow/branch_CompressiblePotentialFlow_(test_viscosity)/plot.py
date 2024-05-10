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

def load_mu_parameters():
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
mu_train, mu_test = load_mu_parameters()
print(" ")
print(len(mu_train), " mu train parameter loaded")
if len(mu_train) > 0: print(mu_train)
print(" ")
print(len(mu_test), " mu test parameter loaded")
if len(mu_test) > 0: print(mu_test)
print(" ")

case_names = ["FOM","ROM","HROM"]
markercolor = ["ob","xr","+g","sc","*m","Dy"]
for j in range(len(mu_train)):
    if os.path.exists("DataBase/Data/" + str(mu_train[j][0]) + ", " + str(mu_train[j][1]) + ".dat"):
        cp_min = 0
        cp_max = 0
        fig = plt.figure()
        fig.set_figwidth(12.0)
        fig.set_figheight(8.0)
        for n, name in enumerate(case_names):
            casename = "Fit"
            if name == "FOM":
                x  = np.loadtxt("DataBase/Data/" + str(mu_train[j][0]) + ", " + str(mu_train[j][1]) + ".dat",usecols=(0,))
                cp = np.loadtxt("DataBase/Data/" + str(mu_train[j][0]) + ", " + str(mu_train[j][1]) + ".dat",usecols=(3,))
                fig = plt.plot(x, cp, markercolor[n], markersize = 3.0, label = name)
                if os.path.exists("DataBase/UncorrectedSolutions/Data/" + str(mu_train[j][0]) + ", " + str(mu_train[j][1]) + ".dat"):
                    aux_x  = np.loadtxt("DataBase/UncorrectedSolutions/Data/" + str(mu_train[j][0]) + ", " + str(mu_train[j][1]) + ".dat",usecols=(0,))
                    aux_cp = np.loadtxt("DataBase/UncorrectedSolutions/Data/" + str(mu_train[j][0]) + ", " + str(mu_train[j][1]) + ".dat",usecols=(3,))
                    fig = plt.plot(aux_x, aux_cp, markercolor[3], markersize = 2.0, label = name + " (uncorrected)")
                if np.min(cp) < cp_min:
                    cp_min = np.min(cp)
                if np.max(cp) > cp_max:
                    cp_max = np.max(cp)

                print("ploting fom case: " + str(mu_train[j][0]) + ", " + str(mu_train[j][1]))

            else:
                if os.path.exists("Data/" + name + "_" + casename + str(mu_train[j][0]) + ", " + str(mu_train[j][1]) + ".dat"):
                    x  = np.loadtxt("Data/" + name + "_" + casename + str(mu_train[j][0]) + ", " + str(mu_train[j][1]) + ".dat",usecols=(0,))
                    cp = np.loadtxt("Data/" + name + "_" + casename + str(mu_train[j][0]) + ", " + str(mu_train[j][1]) + ".dat",usecols=(3,))
                    fig = plt.plot(x, cp, markercolor[n], markersize = 3.0, label = name)
                    if np.min(cp) < cp_min:
                        cp_min = np.min(cp)
                    if np.max(cp) > cp_max:
                        cp_max = np.max(cp)

                    print("ploting ", name)

        fig = plt.title('Cp vs x - ' + casename + " " + str(mu_train[j][0]) + ", " + str(mu_train[j][1]))
        fig = plt.axis([-0.05,1.35,cp_max+0.1,cp_min-0.1])
        fig = plt.ylabel('Cp')
        fig = plt.xlabel('x')
        fig = plt.grid()
        fig = plt.legend()
        fig = plt.tight_layout()
        fig = plt.savefig("Captures/" + casename + str(j) + ".png")
        fig = plt.close('all')

        print(" ")
        print("Train plots ready")
        print(" ")
    else:
        print("The file DataBase/Data/" + str(mu_train[j][0]) + ", " + str(mu_train[j][1]) + ".dat doesn't exist.")


for j in range(len(mu_test)):
    if os.path.exists("DataBase/Data/" + str(mu_test[j][0]) + ", " + str(mu_test[j][1]) + ".dat"):
        cp_min = 0
        cp_max = 0
        fig = plt.figure()
        fig.set_figwidth(12.0)
        fig.set_figheight(8.0)
        for n, name in enumerate(case_names):
            casename = "Test"
            if name == "FOM":
                x  = np.loadtxt("DataBase/Data/" + str(mu_test[j][0]) + ", " + str(mu_test[j][1]) + ".dat",usecols=(0,))
                cp = np.loadtxt("DataBase/Data/" + str(mu_test[j][0]) + ", " + str(mu_test[j][1]) + ".dat",usecols=(3,))
                fig = plt.plot(x, cp, markercolor[n], markersize = 3.0, label = name)
                if np.min(cp) < cp_min:
                    cp_min = np.min(cp)
                if np.max(cp) > cp_max:
                    cp_max = np.max(cp)

                print("ploting fom case:" + str(mu_test[j][0]) + ", " + str(mu_test[j][1]))

            else:
                if os.path.exists("Data/" + name + "_" + casename + str(mu_test[j][0]) + ", " + str(mu_test[j][1]) + ".dat"):
                    x  = np.loadtxt("Data/" + name + "_" + casename + str(mu_test[j][0]) + ", " + str(mu_test[j][1]) + ".dat",usecols=(0,))
                    cp = np.loadtxt("Data/" + name + "_" + casename + str(mu_test[j][0]) + ", " + str(mu_test[j][1]) + ".dat",usecols=(3,))
                    fig = plt.plot(x, cp, markercolor[n], markersize = 3.0, label = name)
                    if np.min(cp) < cp_min:
                        cp_min = np.min(cp)
                    if np.max(cp) > cp_max:
                        cp_max = np.max(cp)

                    print("ploting ", name)

        fig = plt.title('Cp vs x - ' + casename + " " + str(mu_test[j][0]) + ", " + str(mu_test[j][1]))
        fig = plt.axis([-0.05,1.35,cp_max+0.1,cp_min-0.1])
        fig = plt.ylabel('Cp')
        fig = plt.xlabel('x')
        fig = plt.grid()
        fig = plt.legend()
        fig = plt.tight_layout()
        fig = plt.savefig("Captures/" + casename + str(j) + ".png")
        fig = plt.close('all')
    else:
        print("The file DataBase/Data/" + str(mu_test[j][0]) + ", " + str(mu_test[j][1]) + ".dat doesn't exist.")
        
        print(" ")
        print("Test plots ready")
        print(" ")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #