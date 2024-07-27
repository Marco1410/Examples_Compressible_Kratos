import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# load parameters
#
def load_mu_parameters():
    if os.path.exists("mu.npy"):
        mu_values = np.load('mu.npy')
        mu_values =  [mu.tolist() for mu in mu_values]
    return mu_values

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# save parameters
#
def Plot_Cps(mu_list, capture_directory):
    for mu in mu_list:
        case_name = f'{mu[0]}, {mu[1]}'
        capture_filename   = f"{capture_directory}/{case_name}.png"

        #### CP PLOT
        ######################################################################
        cp_min = cp_max = cp_fom = cp_validation = 0
        fig = plt.figure()
        fig.set_figwidth(12.0)
        fig.set_figheight(8.0)

        if capture_directory == 'Validation':
            #### VALIDATION ######
            validation_skin_data_filename = f"reference_data/{case_name}.dat"
            if os.path.exists(validation_skin_data_filename):
                x  = np.loadtxt(validation_skin_data_filename, usecols=(0,))
                cp_validation = np.loadtxt(validation_skin_data_filename, usecols=(1,))
                fig = plt.plot(x, cp_validation, 'r.-', markersize = 5.0, label = 'Validation')
        #### FOM ######
        fom_skin_data_filename = f"FOM_Skin_Data/{case_name}.dat"
        if os.path.exists(fom_skin_data_filename):
            x_fom  = np.loadtxt(fom_skin_data_filename, usecols=(0,))
            cp_fom = np.loadtxt(fom_skin_data_filename, usecols=(3,))
            fig = plt.plot(x_fom, cp_fom, 's', markersize = 5.0, label = 'FOM')
        
        cp_min = np.min([np.min(cp_fom), np.min(cp_validation)])
        cp_max = np.max([np.max(cp_fom), np.max(cp_validation)])
        
        fig = plt.title('Cp vs x')
        fig = plt.axis([-0.05,1.05,cp_max+0.1,cp_min-0.1])
        fig = plt.ylabel('Cp')
        fig = plt.xlabel('x')
        fig = plt.grid()
        fig = plt.legend()
        fig = plt.tight_layout()
        fig = plt.savefig(capture_filename, dpi=400)
        fig = plt.close('all')


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


if __name__ == "__main__":

    mu = load_mu_parameters()

    Plot_Cps(mu, 'Mu_Captures')

    mu_test = []
    mu_test.append([1.0,0.72])
    mu_test.append([1.0,0.73])
    mu_test.append([1.0,0.75])
    mu_test.append([2.0,0.75])

    Plot_Cps(mu_test, 'Validation')