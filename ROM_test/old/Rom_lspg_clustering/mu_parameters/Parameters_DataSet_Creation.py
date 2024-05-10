import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import time as time
import numpy as np
from scipy.stats import qmc

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# get multiple parameters by Halton or LatinHypercube methods
#
def get_multiple_params(number_train_values=0, number_test_values=0, angle=[], mach=[], method='Halton'):
    if np.abs(angle[1]-angle[0])< 1e-2:
        if method == 'Halton':
            sampler = qmc.Halton(d=1)
        elif method == 'LatinHypercube':
            sampler = qmc.LatinHypercube(d=1)
        mu_train = []; mu_test = []
        if number_train_values > 0:
            sample = sampler.random(number_train_values)
            values = qmc.scale(sample, [mach[0]],[mach[1]])
            for i in range(number_train_values):
                #Angle of attack , Mach infinit
                mu_train.append([angle[0], values[i]])
        if number_test_values > 0:
            sample = sampler.random(number_test_values)
            values = qmc.scale(sample, [mach[0]],[mach[1]])
            for i in range(number_test_values):
                #Angle of attack , Mach infinit
                mu_test.append([angle[0], values[i]])
    elif np.abs(mach[1]-mach[0])< 1e-3:
        if method == 'Halton':
            sampler = qmc.Halton(d=1)
        elif method == 'LatinHypercube':
            sampler = qmc.LatinHypercube(d=1)
        mu_train = []; mu_test = []
        if number_train_values > 0:
            sample = sampler.random(number_train_values)
            values = qmc.scale(sample, [angle[0]],[angle[1]])
            for i in range(number_train_values):
                #Angle of attack , Mach infinit
                mu_train.append([values[i], mach[0]])
        if number_test_values > 0:
            sample = sampler.random(number_test_values)
            values = qmc.scale(sample, [angle[0]],[angle[1]])
            for i in range(number_test_values):
                #Angle of attack , Mach infinit
                mu_test.append([values[i], mach[0]])
    else:
        if method == 'Halton':
            sampler = qmc.Halton(d=2)
        elif method == 'LatinHypercube':
            sampler = qmc.LatinHypercube(d=2)
        mu_train = []; mu_test = []
        if number_train_values > 0:
            sample = sampler.random(number_train_values)
            values = qmc.scale(sample, [angle[0],mach[0]], [angle[1],mach[1]])
            for i in range(number_train_values):
                #Angle of attack , Mach infinit
                mu_train.append([values[i,0], values[i,1]])
        if number_test_values > 0:
            sample = sampler.random(number_test_values)
            values = qmc.scale(sample, [angle[0],mach[0]], [angle[1],mach[1]])
            for i in range(number_test_values):
                #Angle of attack , Mach infinit
                mu_test.append([values[i,0], values[i,1]])

    return mu_train, mu_test

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# save parameters
#
def save_mu_parameters(mu_train, name1, mu_test, name2):
    if len(mu_train) > 0:
        np.save(f'{name1}',mu_train)
    if len(mu_test) > 0:
        np.save(f'{name2}',mu_test)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# save parameters
#
def plot_mu_values(mu_train, mu_test):
    if len(mu_train) > 0: plt.plot(mu_train[:,1], mu_train[:,0], 'bs', label="Train Values")
    if len(mu_test) > 0: plt.plot(mu_test[:,1], mu_test[:,0], 'ro', label="Test Values")
    plt.title('Mu Values')
    plt.ylabel('Alpha')
    plt.xlabel('Mach')
    plt.grid(True)
    plt.legend(bbox_to_anchor=(.85, 1.03, 1., .102), loc='upper left', borderaxespad=0.)
    plt.savefig("MuValues.png")
    plt.close('all')


if __name__ == '__main__':

    start_time = time.time()

    NumberofMuTrain = 500
    NumberOfMuTest  = 300

    # Definir rango de valores de mach y angulo de ataque
    mach_range  = [ 0.70, 0.75]
    angle_range = [ 1.00, 2.00]

    mu_train, mu_test = get_multiple_params(number_train_values = NumberofMuTrain,
                        number_test_values  = NumberOfMuTest, 
                        angle  = angle_range, 
                        mach   = mach_range, 
                        method = 'Halton')
    
    mu_train = np.array(mu_train)
    mu_test  = np.array(mu_test)
    
    plot_mu_values(mu_train, mu_test)

    save_mu_parameters(mu_train, 'mu_train', mu_test, 'mu_test')

    exe_time = time.time() - start_time

    print(f' Executing took {round(exe_time, 2)} sec')