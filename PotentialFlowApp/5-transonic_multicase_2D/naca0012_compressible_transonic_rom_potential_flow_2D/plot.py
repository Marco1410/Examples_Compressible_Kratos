import numpy as np
import os
import pickle
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

archivo = open('Data/mu_train.dat', 'rb')
lista = pickle.load(archivo)
mu_train = np.asarray(lista)
archivo.close()
archivo = open('Data/mu_test.dat', 'rb')
lista = pickle.load(archivo)
mu_test = np.asarray(lista)
archivo.close()
   
case_names = ["FOM","ROM","HROM"]
markercolor = ["ob","xr","+g"]
database_path = "../TransonicDataBase/"
for j in range(len(mu_train)):
    cp_min = 0
    cp_max = 0
    fig = plt.figure()
    fig.set_figwidth(12.0)
    fig.set_figheight(8.0)
    for n, name in enumerate(case_names):
        if name == "FOM":
            x  = np.loadtxt(database_path+"Data/"+name+"_"+str(mu_train[j][0])+"_"+str(mu_train[j][1])+".dat",usecols=(0,))
            cp = np.loadtxt(database_path+"Data/"+name+"_"+str(mu_train[j][0])+"_"+str(mu_train[j][1])+".dat",usecols=(3,))
            fig = plt.plot(x, cp, markercolor[n], markersize = 3.0, label = name)
            if np.min(cp) < cp_min:
                cp_min = np.min(cp)
            if np.max(cp) > cp_max:
                cp_max = np.max(cp)
        elif os.path.exists("Data/"+name+"_"+"train"+"_"+str(j)+".dat"):
            x  = np.loadtxt("Data/"+name+"_"+"train"+"_"+str(j)+".dat",usecols=(0,))
            cp = np.loadtxt("Data/"+name+"_"+"train"+"_"+str(j)+".dat",usecols=(3,))
            fig = plt.plot(x, cp, markercolor[n], markersize = 3.0, label = name)
            if np.min(cp) < cp_min:
                cp_min = np.min(cp)
            if np.max(cp) > cp_max:
                cp_max = np.max(cp)
    fig = plt.title('Cp vs x - ' + "angle: " + str(mu_train[j][0]) + "º " + "mach: " + str(mu_train[j][1]))
    fig = plt.axis([-0.05,1.35,cp_max+0.1,cp_min-0.1])
    fig = plt.ylabel('Cp')
    fig = plt.xlabel('x')
    fig = plt.grid()
    fig = plt.legend()
    fig = plt.tight_layout()
    fig = plt.savefig("Captures/Simulation_train_" + str(j) + "_A_" + str(mu_train[j][0]) + "_M_" + str(mu_train[j][1])+".png")
    fig = plt.close('all')

for j in range(len(mu_test)):
    cp_min = 0
    cp_max = 0
    fig = plt.figure()
    fig.set_figwidth(12.0)
    fig.set_figheight(8.0)
    for n, name in enumerate(case_names):
        if name == "FOM":
            x  = np.loadtxt(database_path+"Data/"+name+"_"+str(mu_test[j][0])+"_"+str(mu_test[j][1])+".dat",usecols=(0,))
            cp = np.loadtxt(database_path+"Data/"+name+"_"+str(mu_test[j][0])+"_"+str(mu_test[j][1])+".dat",usecols=(3,))
            fig = plt.plot(x, cp, markercolor[n], markersize = 3.0, label = name)
            if np.min(cp) < cp_min:
                cp_min = np.min(cp)
            if np.max(cp) > cp_max:
                cp_max = np.max(cp)
        elif os.path.exists("Data/"+name+"_"+"test"+"_"+str(j)+".dat"):
            x  = np.loadtxt("Data/"+name+"_"+"test"+"_"+str(j)+".dat",usecols=(0,))
            cp = np.loadtxt("Data/"+name+"_"+"test"+"_"+str(j)+".dat",usecols=(3,))
            fig = plt.plot(x, cp, markercolor[n], markersize = 3.0, label = name)
            if np.min(cp) < cp_min:
                cp_min = np.min(cp)
            if np.max(cp) > cp_max:
                cp_max = np.max(cp)
    fig = plt.title('Cp vs x - ' + "angle: " + str(mu_test[j][0]) + "º " + "mach: " + str(mu_test[j][1]))
    fig = plt.axis([-0.05,1.35,cp_max+0.1,cp_min-0.1])
    fig = plt.ylabel('Cp')
    fig = plt.xlabel('x')
    fig = plt.grid()
    fig = plt.legend()
    fig = plt.tight_layout()
    fig = plt.savefig("Captures/Simulation_test_" + str(j) + "_A_" + str(mu_test[j][0]) + "_M_" + str(mu_test[j][1])+".png")
    fig = plt.close('all')