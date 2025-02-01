from ezyrb import ReducedOrderModel as ROM
from ezyrb import RBF, POD, Database
import numpy as np
import os




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# load parameters
#
def load_mu_parameters(name):
    filename = f'{name}.npy'
    if os.path.exists(filename):
        mu_npy = np.load(filename)
        mu =  [mu.tolist() for mu in mu_npy]
    else:
        mu = []
    return mu

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# RBF prediction and error
#
def RBF_prediction( mu_train                 = [None],
                    mu_test                  = [None],
                    mu_train_not_scaled      = [None],
                    mu_test_not_scaled       = [None],
                    full                     = False):

    #### CLUSTERING DATA
    #################################################################################################
    parameters = np.array(mu_train_not_scaled)

    snapshots = []
    for mu in mu_train:
        file = f'{mu[0]}, {mu[1]}'
        if full:
            snapshots.append(np.array(np.load(f'FOM_Snapshots/{file}.npy')).reshape(-1,1))
        else:
            snapshots.append(np.array(np.loadtxt(f'FOM_Skin_Data/{file}.dat', usecols=(3,))).reshape(-1,1))

    snapshots = np.block(snapshots)

    #### RBF TRAINING
    #################################################################################################
    db = Database(parameters, snapshots.T)
    pod = POD()
    rbf = RBF()
    rom = ROM(db, pod, rbf).fit()

    if len(mu_test) > 0:
        #### PREDICTION OF TEST
        #################################################################################################
        interpolated_solutions_list = [rom.predict([element]).snapshots_matrix for element in mu_test_not_scaled]

        for i, solution in enumerate(interpolated_solutions_list):
            if full:
                np.save(f"RBF_Snapshots/{mu_test[i][0]}, {mu_test[i][1]}", solution.T)
            else:
                np.save(f"RBF_Skin_Data/{mu_test[i][0]}, {mu_test[i][1]}", solution.T)

def RBF_error_estimation(mu_train, mu_test, full=False):
    if len(mu_train) > 0:
        approximation_error = 0.0
        FOM_model = []; RBF_model = []
        for mu in mu_train:
            if full:
                FOM_model.append(np.array(np.load(f'FOM_Snapshots/{mu[0]}, {mu[1]}.npy')).reshape(-1,1))
                RBF_model.append(np.array(np.load(f"RBF_Snapshots/{mu[0]}, {mu[1]}.npy")).reshape(-1,1))
            else:
                FOM_model.append(np.array(np.loadtxt(f'FOM_Skin_Data/{mu[0]}, {mu[1]}.dat', usecols=(3,))).reshape(-1,1))
                RBF_model.append(np.array(np.load(f"RBF_Skin_Data/{mu[0]}, {mu[1]}.npy")).reshape(-1,1))
        FOM_model = np.block(FOM_model)
        RBF_model = np.block(RBF_model)
        training_approximation_error = np.linalg.norm(FOM_model - RBF_model)/np.linalg.norm(FOM_model)
        print(f'RBF training approximation error: {training_approximation_error:.2E}')

    if len(mu_test)>0:
        if full:
            path = f"RBF_Snapshots/{mu_test[0][0]}, {mu_test[0][1]}.npy"
        else:
            path = f'RBF_Skin_Data/{mu_test[0][0]}, {mu_test[0][1]}.npy'

        if os.path.exists(path):
            approximation_error = 0.0
            FOM_model = []; RBF_model_interpolation = []
            for mu in mu_test:
                if full:
                    FOM_model.append(np.array(np.load(f'FOM_Snapshots/{mu[0]}, {mu[1]}.npy')).reshape(-1,1))
                    RBF_model_interpolation.append(np.array(np.load(f"RBF_Snapshots/{mu[0]}, {mu[1]}.npy")).reshape(-1,1))
                else:
                    FOM_model.append(np.array(np.loadtxt(f'FOM_Skin_Data/{mu[0]}, {mu[1]}.dat', usecols=(3,))).reshape(-1,1))
                    RBF_model_interpolation.append(np.array(np.load(f"RBF_Skin_Data/{mu[0]}, {mu[1]}.npy")).reshape(-1,1))
            FOM_model = np.block(FOM_model)
            RBF_model_interpolation = np.block(RBF_model_interpolation)
            approximation_error = np.linalg.norm(FOM_model - RBF_model_interpolation)/np.linalg.norm(FOM_model)
            print(f'RBF interpolation approximation error: {approximation_error:.2E}')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



if __name__ == "__main__":

    full = True

    strategies = ['galerkin']
    n = 0

    mu_train      = load_mu_parameters(f'Mu_history/{n}_{strategies[0]}_mu_train')
    mu_train_not_scaled = load_mu_parameters(f'Mu_history/{n}_{strategies[0]}_mu_train_not_scaled')

    mu_test       = load_mu_parameters(f'Mu_history/{n}_{strategies[0]}_mu_test')
    mu_test_not_scaled = load_mu_parameters(f'Mu_history/{n}_{strategies[0]}_mu_test_not_scaled')

    mu_validation = load_mu_parameters(f'Mu_history/mu_validation')
    mu_validation_not_scaled = load_mu_parameters(f'Mu_history/mu_validation_not_scaled')

    if len(mu_train) >= 3:
        RBF_prediction(mu_train = mu_train, mu_train_not_scaled = mu_train_not_scaled, 
                    mu_test = mu_train + mu_test, mu_test_not_scaled  = mu_train_not_scaled + mu_test_not_scaled, full=full)
        RBF_prediction(mu_train = mu_train, mu_train_not_scaled = mu_train_not_scaled, 
                    mu_test = mu_validation, mu_test_not_scaled  = mu_validation_not_scaled, full=full)


    if len(mu_train) >= 3:
        print('::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
        RBF_error_estimation(mu_train, mu_test, full=full)
        print('::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
        print('Validation Error::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
        print('::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
        RBF_error_estimation(mu_train, mu_validation, full=full)
        print('::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::')