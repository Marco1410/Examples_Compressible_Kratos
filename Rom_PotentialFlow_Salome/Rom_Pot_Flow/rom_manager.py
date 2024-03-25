import os
import random
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import pickle
import numpy as np
from scipy.stats import qmc
import time as time
import KratosMultiphysics
import KratosMultiphysics.kratos_utilities
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
from local_rom_manager import LocalRomManager

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def CustomizeSimulation(cls, global_model, parameters):

    class CustomSimulation(cls):

        def __init__(self, model,project_parameters, custom_param = None):
            super().__init__(model,project_parameters)
            self.custom_param  = custom_param

            # if self._GetSimulationName() == "::[ROM Simulation]:: ":
                #just in case the residual not converge
                # parameters["solver_settings"]["solving_strategy_settings"]["type"].SetString("newton_raphson")
                # parameters["solver_settings"]["maximum_iterations"].SetInt(150)
                # parameters["solver_settings"]["relative_tolerance"].SetDouble(1e-12)
                # parameters["solver_settings"]["scheme_settings"]["update_transonic_tolerance"].SetDouble(1e-2)
        
        def Initialize(self):
            super().Initialize()

        def InitializeSolutionStep(self):
            super().InitializeSolutionStep()

            set_initial_conditions_rom_from_fom = False
            # To set initial conditions:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            if self._GetSimulationName() == "::[ROM Simulation]:: ":
                
                self._GetSolver()._GetBuilderAndSolver().SetCorrectWithFOM(True)

                if parameters["output_processes"].Has("vtk_output") and set_initial_conditions_rom_from_fom:
                    nametype = parameters["output_processes"]["vtk_output"][0]["Parameters"]["output_path"].GetString()
                    simulation_name = nametype.split('/')[1].split('_')[1].removeprefix("Fit").removeprefix("Test")
                    data_name = "DataBase/" + simulation_name + ".npy"
                    data = np.load(data_name)
                    node_ids = np.load("rom_data/NodeIds.npy")

                    for node in self.model["MainModelPart"].Nodes:
                        offset = np.where(node_ids == node.Id)[0][0]*2
                        node.SetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL, data[offset+1,0])
                        node.SetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL, data[offset,0])

                # data_name = "initial_solution_cluster_1.npy"
                # u = np.load(data_name)
                # phi  = np.load("rom_data/RightBasisMatrix.npy")
                # data = phi@(phi.T@u)
                # node_ids = np.load("rom_data/NodeIds.npy")

                # for node in self.model["MainModelPart"].Nodes:
                #     offset = np.where(node_ids == node.Id)[0][0]*2
                #     node.SetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL, data[offset+1])
                #     node.SetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL, data[offset])

        def FinalizeSolutionStep(self):
            super().FinalizeSolutionStep()

            if parameters["output_processes"].Has("vtk_output"):
                nametype = parameters["output_processes"]["vtk_output"][0]["Parameters"]["output_path"].GetString()

                if self._GetSimulationName() == "Analysis":
                    simulation_name = nametype.removeprefix("DataBase/Results/").removeprefix("Results/FOM_Test")
                    skin_data = "DataBase/Data/" + simulation_name + ".dat"
                    fout = open(skin_data,'w')
                else:
                    simulation_name = nametype.split('/')[1]
                    skin_data = "Data/" + simulation_name + ".dat"
                    fout = open(skin_data,'w')

                modelpart = self.model["MainModelPart.Body2D_Body"]
                for node in modelpart.Nodes:
                    x = node.X ; y = node.Y ; z = node.Z
                    cp = node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
                    fout.write("%s %s %s %s\n" %(x,y,z,cp))
                fout.close()

            if self._GetSimulationName() == "::[ROM Simulation]:: ":
                if parameters["output_processes"].Has("vtk_output"):
                    nametype = parameters["output_processes"]["vtk_output"][0]["Parameters"]["output_path"].GetString()
                    simulation_name = nametype.split('/')[1].split('_')[1].removeprefix("Fit").removeprefix("Test")
                    if not os.path.exists("DataBase/aux_data"):
                        os.mkdir("DataBase/aux_data")
                    aux_data_name = "DataBase/aux_data/" + simulation_name + ".npy"
                    if np.array(self._GetSolver()._GetBuilderAndSolver().GetIntermediateSolutionsMatrix()).ndim >=1: 
                        if not os.path.exists(aux_data_name):
                            np.save(aux_data_name, np.array(self._GetSolver()._GetBuilderAndSolver().GetIntermediateSolutionsMatrix()))
                        else:
                            aux_list = np.concatenate((np.load(aux_data_name),np.array(self._GetSolver()._GetBuilderAndSolver().GetIntermediateSolutionsMatrix())), axis=1)
                            np.save(aux_data_name, aux_list)

    return CustomSimulation(global_model, parameters)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def UpdateProjectParameters(parameters, mu=None):

    angle_of_attack        = mu[0]
    mach_infinity          = mu[1]

    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(angle_of_attack)
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach_infinity)

    return parameters

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# get multiple parameters
#
def get_multiple_params_by_Halton_sequence(number_train_values,number_test_values,angle,mach,fix_corners_of_parametric_space):
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_train.dat')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_test.dat')
    if np.abs(angle[1]-angle[0])< 1e-2:
        if fix_corners_of_parametric_space and number_train_values < 2:
            print("Setting number of values to 2.")
            number_train_values = 2
        sampler = qmc.Halton(d=1)
        # sampler = qmc.LatinHypercube(d=1)
        mu_train = []; mu_test = []
        if number_train_values > 0:
            sample = sampler.random(number_train_values)
            values = qmc.scale(sample, [mach[0]],[mach[1]])
            if fix_corners_of_parametric_space and number_train_values >= 2:
                values[0] = mach[0]
                values[number_train_values-1] = mach[1]
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
        if fix_corners_of_parametric_space and number_train_values < 2:
            print("Setting number of values to 2.")
            number_train_values = 2
        sampler = qmc.Halton(d=1)
        # sampler = qmc.LatinHypercube(d=1)
        mu_train = []; mu_test = []
        if number_train_values > 0:
            sample = sampler.random(number_train_values)
            values = qmc.scale(sample, [angle[0]],[angle[1]])
            if fix_corners_of_parametric_space and number_train_values >= 2:
                values[0] = angle[0]
                values[number_train_values-1] = angle[1]
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
        if fix_corners_of_parametric_space and number_train_values < 4:
            print("Setting number of values to 4.")
            number_train_values = 4
        sampler = qmc.Halton(d=2)
        # sampler = qmc.LatinHypercube(d=2)
        mu_train = []; mu_test = []
        if number_train_values > 0:
            sample = sampler.random(number_train_values)
            values = qmc.scale(sample, [angle[0],mach[0]], [angle[1],mach[1]])
            if fix_corners_of_parametric_space and number_train_values >= 4:
                values[0,0] = angle[0]
                values[0,1] = mach[0]
                values[1,0] = angle[0]
                values[1,1] = mach[1]
                values[number_train_values-1,0] = angle[1]
                values[number_train_values-1,1] = mach[1]
                values[number_train_values-2,0] = angle[1]
                values[number_train_values-2,1] = mach[0]
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

def plot_mu_values(mu_train,mu_test):
    mu_train_a = np.zeros(len(mu_train))
    mu_train_m = np.zeros(len(mu_train))
    mu_test_a  = np.zeros(len(mu_test))
    mu_test_m  = np.zeros(len(mu_test))
    for i in range(len(mu_train)):
        mu_train_a[i] = mu_train[i][0]
        mu_train_m[i] = mu_train[i][1]
    for i in range(len(mu_test)):
        mu_test_a[i] = mu_test[i][0]
        mu_test_m[i] = mu_test[i][1]
    plt.plot(mu_train_m, mu_train_a, 'bs', label="Train Values")
    plt.plot(mu_test_m, mu_test_a, 'ro', label="Test Values")
    plt.title('Mu Values')
    plt.ylabel('Alpha')
    plt.xlabel('Mach')
    plt.grid(True)
    plt.legend(bbox_to_anchor=(.85, 1.03, 1., .102), loc='upper left', borderaxespad=0.)
    plt.savefig("MuValues.png")
    plt.close('all')

def plot_mu_values_with_errors():
    mu_train, mu_test = load_mu_parameters(with_errors=True)
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
#
# plot
#

def plot_Cps():
    mu_train, mu_test = load_mu_parameters()
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
                    fig = plt.plot(x, cp, markercolor[n], markersize = 2.0, label = name)
                    if np.min(cp) < cp_min:
                        cp_min = np.min(cp)
                    if np.max(cp) > cp_max:
                        cp_max = np.max(cp)
                else:
                    if os.path.exists("Data/" + name + "_" + casename + str(mu_train[j][0]) + ", " + str(mu_train[j][1]) + ".dat"):
                        x  = np.loadtxt("Data/" + name + "_" + casename + str(mu_train[j][0]) + ", " + str(mu_train[j][1]) + ".dat",usecols=(0,))
                        cp = np.loadtxt("Data/" + name + "_" + casename + str(mu_train[j][0]) + ", " + str(mu_train[j][1]) + ".dat",usecols=(3,))
                        fig = plt.plot(x, cp, markercolor[n], markersize = 2.0, label = name)
                        if np.min(cp) < cp_min:
                            cp_min = np.min(cp)
                        if np.max(cp) > cp_max:
                            cp_max = np.max(cp)
            fig = plt.title('Cp vs x - ' + casename + " " + str(mu_train[j][0]) + ", " + str(mu_train[j][1]))
            fig = plt.axis([-0.05,1.35,cp_max+0.1,cp_min-0.1])
            fig = plt.ylabel('Cp')
            fig = plt.xlabel('x')
            fig = plt.grid()
            fig = plt.legend()
            fig = plt.tight_layout()
            fig = plt.savefig("Captures/" + casename + str(j) + ".png")
            fig = plt.close('all')
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
                    fig = plt.plot(x, cp, markercolor[n], markersize = 2.0, label = name)
                    if np.min(cp) < cp_min:
                        cp_min = np.min(cp)
                    if np.max(cp) > cp_max:
                        cp_max = np.max(cp)
                else:
                    if os.path.exists("Data/" + name + "_" + casename + str(mu_test[j][0]) + ", " + str(mu_test[j][1]) + ".dat"):
                        x  = np.loadtxt("Data/" + name + "_" + casename + str(mu_test[j][0]) + ", " + str(mu_test[j][1]) + ".dat",usecols=(0,))
                        cp = np.loadtxt("Data/" + name + "_" + casename + str(mu_test[j][0]) + ", " + str(mu_test[j][1]) + ".dat",usecols=(3,))
                        fig = plt.plot(x, cp, markercolor[n], markersize = 2.0, label = name)
                        if np.min(cp) < cp_min:
                            cp_min = np.min(cp)
                        if np.max(cp) > cp_max:
                            cp_max = np.max(cp)
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

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# save / load parameters
#
def save_mu_parameters(mu_train, mu_test, name1, name2):
    if len(mu_train) > 0:
        archivo = open(name1+'.dat', 'wb')
        pickle.dump(mu_train, archivo)
        archivo.close()
    if len(mu_test) > 0:
        archivo = open(name2+'.dat', 'wb')
        pickle.dump(mu_test, archivo)
        archivo.close()

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

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
def CleanFolder(database=False):
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Results')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Data')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Captures')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('MuValues.png')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('MuValuesWithErrors.png')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_train_errors.dat')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_test_errors.dat')
    os.mkdir("Results")
    os.mkdir("Data")
    os.mkdir("Captures")
    if database:
        KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('DataBase')

def CleanToTest(database=False):
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Results')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Data')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Captures')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('MuValues.png')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('MuValuesWithErrors.png')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_train_errors.dat')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_test_errors.dat')
    os.mkdir("Results")
    os.mkdir("Data")
    os.mkdir("Captures")
    if database:
        KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('DataBase')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def GetRomManagerParameters():
    general_rom_manager_parameters = KratosMultiphysics.Parameters("""{
            "rom_stages_to_train" : ["ROM"],      // ["ROM","HROM"]
            "rom_stages_to_test"  : ["ROM"],      // ["ROM","HROM"]
            "paralellism" : null,                        // null, TODO: add "compss"
            "projection_strategy": "galerkin",           // "lspg", "galerkin", "petrov_galerkin"
            "assembling_strategy": "global",             // "global", "elemental"
            "save_gid_output": false,                     // false, true #if true, it must exits previously in the ProjectParameters.json
            "save_vtk_output": true,
            "output_name": "mu",                         // "id" , "mu"
            "ROM":{
                "svd_truncation_tolerance": 1e-12,
                "model_part_name": "MainModelPart",
                "nodal_unknowns": ["VELOCITY_POTENTIAL","AUXILIARY_VELOCITY_POTENTIAL"], // Main unknowns. Snapshots are taken from these
                "rom_basis_output_format": "numpy",                                       // "json" "numpy"
                "rom_basis_output_name": "RomParameters",
                "rom_basis_output_folder": "rom_data",
                "snapshots_control_type": "step",                                        // "step", "time"
                "snapshots_interval": 1,
                "galerkin_rom_bns_settings": {
                    "monotonicity_preserving": false
                },
                "lspg_rom_bns_settings": {
                    "train_petrov_galerkin": false,
                    "basis_strategy": "reactions",                        // 'residuals', 'jacobian', 'reactions'
                    "include_phi": false,
                    "svd_truncation_tolerance": 1e-12,
                    "solving_technique": "normal_equations",              // 'normal_equations', 'qr_decomposition'
                    "monotonicity_preserving": false
                },
                "petrov_galerkin_rom_bns_settings": {
                    "monotonicity_preserving": false
                }
            },
            "HROM":{
                "element_selection_type": "empirical_cubature",
                "initial_candidate_elements_model_part_list": [],
                "element_selection_svd_truncation_tolerance": 1e-12,
                "create_hrom_visualization_model_part" : false,
                "echo_level" : 0
            }
        }""")
    return general_rom_manager_parameters

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def create_salome_mesh(salome_script_name):
    import subprocess

    start_time = time.time()

    salome_cmd = "salome -t python"
    salome_exe = " ".join([salome_cmd, salome_script_name])
    sp = subprocess.Popen(["/bin/bash", "-i", "-c", salome_exe])
    sp.communicate()

    exe_time = time.time() - start_time

    print('Executing SALOME with "' + salome_script_name +
        '" took ' + str(round(exe_time, 2)) + ' sec')
    
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

if __name__ == "__main__":

    # if not os.path.exists("salome_files"):
    #     create_salome_mesh("generate_naca0012_2d.py")

    start_time = time.time()

    NumberofMuTrain = 25
    NumberOfMuTest  = 0

    mu_from_database                   = False
    load_old_mu_parameters             = True
    fix_corners_train_parametric_space = False
    only_test                          = False
    update_database                    = False
    use_non_linear_steps               = True
    use_aux_step                       = True

    #CLUSTERING SETTINGS
    clustering_method                  = "k-means" # PEBL - k-means - PEBL_full
    n_clusters                         = 8
    bisection_tolerance                = 0.15
    POD_tolerance                      = 1e-12

    dir = "DataBase"

    # Definir rango de valores de mach y angulo de ataque
    mach_range  = [ 0.65, 0.75]
    angle_range = [ 0.00, 2.50]

    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    if only_test:
        CleanToTest(update_database)
    else:
        CleanFolder(update_database)

    if mu_from_database:
        data = []
        with os.scandir(dir) as files:
            files = [file.name for file in files if file.is_file() and file.name.endswith('.npy')]
        for file in files:
            name = file.removesuffix(".npy")
            angle = np.double(name.split(",")[0])
            mach  = np.double(name.split(",")[1].removeprefix(" "))

            if (    angle <= angle_range[1]
                and angle >= angle_range[0] 
                and mach  <=  mach_range[1]
                and mach  >=  mach_range[0]):
                data.append([angle, mach])

        if len(data) >= NumberofMuTrain:
            mu_train = random.sample(data, NumberofMuTrain)
        else:
            print("mu_train size: ", len(data))
            input("pause PRESS to continue")
            mu_train = data
        if len(data) >= NumberOfMuTest:
            mu_test  = random.sample(data, NumberOfMuTest)
        else:
            print("mu_test size: ", len(data))
            input("pause PRESS to continue")
            mu_test = data
    else:
        if load_old_mu_parameters:
            mu_train, mu_test = load_mu_parameters()
        else:
            mu_train, mu_test = get_multiple_params_by_Halton_sequence(NumberofMuTrain, NumberOfMuTest, angle_range, mach_range,
                                                            fix_corners_train_parametric_space)
  
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    save_mu_parameters(mu_train,mu_test,"mu_train","mu_test")
    plot_mu_values(mu_train,mu_test)
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    general_rom_manager_parameters = GetRomManagerParameters()
    project_parameters_name = "ProjectParameters.json"
    rom_manager = LocalRomManager(project_parameters_name,general_rom_manager_parameters,CustomizeSimulation,UpdateProjectParameters)
    mu_train_errors = rom_manager.Fit(mu_train,use_non_linear_steps, clustering_method, n_clusters, bisection_tolerance, POD_tolerance, use_aux_step)
    mu_test_errors  = rom_manager.Test(mu_test)
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    save_mu_parameters(mu_train_errors,mu_test_errors,"mu_train_errors","mu_test_errors")
    plot_mu_values_with_errors()
    rom_manager.PrintErrors()
    plot_Cps()

    exe_time = time.time() - start_time

    print('Executing ROM Manager with "rom_manager.py" took ' + str(round(exe_time, 2)) + ' sec')