import os
import random
import importlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
from ezyrb import Database
from ezyrb import RBF, POD
from ezyrb import ReducedOrderModel as ROM
import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp 
from KratosMultiphysics.RomApplication.calculate_rom_basis_output_process import CalculateRomBasisOutputProcess


def load_mu_parameters():
    if os.path.exists("mu_train.npy") and os.path.exists("mu_test.npy"):
        archivo = 'mu_train.npy'
        mu_train = np.load(archivo)
        archivo = 'mu_test.npy'
        mu_test = np.load(archivo)
    elif os.path.exists("mu_train.npy"):
        archivo = 'mu_train.npy'
        mu_train = np.load(archivo)
        mu_test = []
    elif os.path.exists("mu_test.npy"):
        archivo = 'mu_test.npy'
        mu_test = np.load(archivo)
        mu_train = []
    return list(mu_train), list(mu_test)

def error_estimation(mu_train, mu_list):
    approximation_error = 0.0
    FOM_model = []; RBF_model = []
    # LaunchSimulations(mu_list, False)
    for mu in mu_train:
        FOM_model.append(np.load(f'DataBase/{mu[0]}, {mu[1]}.npy'))
        RBF_model.append(np.load(f"Interpolated_data/{mu[0]}, {mu[1]}.npy").T)
    FOM_model = np.block(FOM_model)
    RBF_model = np.block(RBF_model)
    training_approximation_error = np.linalg.norm(FOM_model - RBF_model)/np.linalg.norm(FOM_model)
    print(f'RBF training approximation error: {training_approximation_error*100} %')

    approximation_error = 0.0
    FOM_model = []; RBF_model = []; RBF_model_interpolation = []
    for mu in mu_list:
        FOM_model.append(np.load(f'DataBase/{mu[0]}, {mu[1]}.npy'))
        RBF_model_interpolation.append(np.load(f"Interpolated_data/{mu[0]}, {mu[1]}.npy").T)
        # RBF_model.append(np.load(f"NR/{mu[0]}, {mu[1]}.npy"))
    FOM_model = np.block(FOM_model)
    # RBF_model = np.block(RBF_model)
    RBF_model_interpolation = np.block(RBF_model_interpolation)
    # approximation_error = np.linalg.norm(FOM_model - RBF_model)/np.linalg.norm(FOM_model)
    # print(f'Interpolation approximation error NR: {approximation_error*100} %')
    approximation_error = np.linalg.norm(FOM_model - RBF_model_interpolation)/np.linalg.norm(FOM_model)
    print(f'Interpolation approximation error: {approximation_error*100} %')
    return approximation_error*100

def LaunchSimulations(mu_test, test):
    if not os.path.exists("DataBase"):
        os.mkdir("DataBase")
    if not os.path.exists("NR"):
        os.mkdir("NR")
    with open('ProjectParameters.json','r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    for mu in mu_test:
        if not os.path.exists(f'DataBase/{mu[0]}, {mu[1]}.npy') and test:
            model = KratosMultiphysics.Model()
            parameters_copy = UpdateProjectParameters(parameters.Clone(), mu)
            analysis_stage_class = _GetAnalysisStageClass(parameters_copy)
            simulation = CustomizeSimulation(analysis_stage_class,model,parameters_copy, test)
            simulation.Run()
            for process in simulation._GetListOfOutputProcesses():
                if isinstance(process, CalculateRomBasisOutputProcess):
                    BasisOutputProcess = process
            np.save(f'DataBase/{mu[0]}, {mu[1]}', BasisOutputProcess._GetSnapshotsMatrix())

        elif not os.path.exists(f'NR/{mu[0]}, {mu[1]}.npy') and not test:
            model = KratosMultiphysics.Model()
            parameters_copy = UpdateProjectParameters(parameters.Clone(), mu)
            analysis_stage_class = _GetAnalysisStageClass(parameters_copy)
            simulation = CustomizeSimulation(analysis_stage_class,model,parameters_copy, test)
            simulation.Run()
            for process in simulation._GetListOfOutputProcesses():
                if isinstance(process, CalculateRomBasisOutputProcess):
                    BasisOutputProcess = process
            np.save(f'NR/{mu[0]}, {mu[1]}', BasisOutputProcess._GetSnapshotsMatrix())

def CustomizeSimulation(cls, global_model, parameters, test=False):
    class CustomSimulation(cls):
        def __init__(self, model,project_parameters, custom_param = None):
            super().__init__(model,project_parameters)
            self.custom_param  = custom_param
        def Initialize(self):
            super().Initialize()
        def InitializeSolutionStep(self):
            super().InitializeSolutionStep()
            
            if not test:
                angle = parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].GetDouble()
                mach  = parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].GetDouble()
                data_name = f"Interpolated_data/{angle}, {mach}.npy"
                u = np.load(data_name).T

                model_part = self.model["MainModelPart"]
                for node in model_part.Nodes:
                    offset = np.where(np.arange(1,model_part.NumberOfNodes()+1, dtype=int) == node.Id)[0][0]*2
                    node.SetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL, u[offset+1])
                    node.SetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL, u[offset])
            
        def FinalizeSolutionStep(self):
            super().FinalizeSolutionStep()

    return CustomSimulation(global_model, parameters)

def UpdateProjectParameters(parameters, mu=None):
    angle_of_attack        = mu[0]
    mach_infinity          = mu[1]
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(angle_of_attack)
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach_infinity)
    return parameters

def _GetAnalysisStageClass(parameters):
    analysis_stage_module_name = parameters["analysis_stage"].GetString()
    analysis_stage_class_name = analysis_stage_module_name.split('.')[-1]
    analysis_stage_class_name = ''.join(x.title() for x in analysis_stage_class_name.split('_'))
    analysis_stage_module = importlib.import_module(analysis_stage_module_name)
    analysis_stage_class = getattr(analysis_stage_module, analysis_stage_class_name)
    return analysis_stage_class

if __name__ == '__main__':

    mu_train, mu_test = load_mu_parameters()

    # NumberOfMuTrain = 10
    # mu_train  = random.sample(mu_train, NumberOfMuTrain)

    # NumberOfMuTest  = 10
    # mu_test   = random.sample(mu_test, NumberOfMuTest)

    # with os.scandir(f'DataBase') as files:
    #     files = [file.name.removesuffix('.npy') for file in files if file.is_file() and file.name.endswith('.npy')]

    # mu_train = [[np.double(file.split(',')[0]), np.double(file.split(',')[1].removeprefix(' '))] for file in files]
    # mu_train = [np.array(mu) for mu in np.array(mu_train) if mu not in np.array(mu_test)]

    test_error_estimation = 1

    tolerance = 1e-5

    while test_error_estimation > tolerance: 

        #### DATA
        #################################################################################################
        snapshots = []
        for mu in list(mu_train):
            file = f'{mu[0]}, {mu[1]}.npy'
            snapshots.append(np.load(f'DataBase/{file}'))
        snapshots = np.block(snapshots)
        np.save(f'SnapshotsMatrix_conv', snapshots)

        print(f'Snapshots Matrix shape: {snapshots.shape}')
        print(f'Parameters shape: {np.array(mu_train).shape}')

        #### RBF TRAINING
        #################################################################################################
        db = Database(np.array(mu_train), snapshots.T)
        pod = POD()
        rbf = RBF()
        rom = ROM(db, pod, rbf).fit()

        # for pt, error in zip(parameters, rom.loo_error()):
        #     print(pt, error)

        training_solutions_list = [rom.predict([element]).snapshots_matrix for element in mu_train]
        if not os.path.exists("Interpolated_data"):
            os.mkdir("Interpolated_data")
        for i, solution in enumerate(training_solutions_list):
            np.save(f"Interpolated_data/{mu_train[i][0]}, {mu_train[i][1]}.npy", solution)

        solutions_list = [rom.predict([element]).snapshots_matrix for element in mu_test]
        if not os.path.exists("Interpolated_data"):
            os.mkdir("Interpolated_data")
        for i, solution in enumerate(solutions_list):
            np.save(f"Interpolated_data/{mu_test[i][0]}, {mu_test[i][1]}.npy", solution)
        LaunchSimulations(mu_test, True)

        test_error_estimation = error_estimation(mu_train, mu_test)
        
        if test_error_estimation > tolerance:
            new_mu = rom.optimal_mu()

            print(f'New mu: {new_mu}')
            mu_train.append([new_mu[0][0], new_mu[0][1]])

            if not os.path.exists(f'DataBase/{new_mu[0][0]}, {new_mu[0][1]}.npy'):
                LaunchSimulations(new_mu, True)

            np.save('new_mu_train', mu_train)
