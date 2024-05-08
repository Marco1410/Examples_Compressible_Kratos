import os
import random
import importlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
from ezyrb import Database
from ezyrb import RBF, POD
from ezyrb import ReducedOrderModel as ROM
import torch
from ezyrb import Linear, ANN, Snapshot, Parameter
from ezyrb.plugin import AutomaticShiftSnapshots
import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp 
import KratosMultiphysics.MedApplication as KratosMED
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
    return np.array(mu_train), np.array(mu_test)

def error_estimation(mu_train, mu_list):
    approximation_error = 0.0
    FOM_model = []; RBF_model = []
    for mu in mu_train:
        FOM_model.append(np.load(f'DataBase/{mu[0]}, {mu[1]}.npy'))
        RBF_model.append(np.load(f"Interpolated_data/{mu[0]}, {mu[1]}.npy").T)
    FOM_model = np.block(FOM_model)
    RBF_model = np.block(RBF_model)
    training_approximation_error = np.linalg.norm(FOM_model - RBF_model)/np.linalg.norm(FOM_model)
    print(f'RBF training approximation error: {training_approximation_error*100} %')

    approximation_error = 0.0
    FOM_model = []; RBF_model = []
    for mu in mu_list:
        FOM_model.append(np.load(f'DataBase/{mu[0]}, {mu[1]}.npy'))
        RBF_model.append(np.load(f"Interpolated_data/{mu[0]}, {mu[1]}.npy").T)
    FOM_model = np.block(FOM_model)
    RBF_model = np.block(RBF_model)
    approximation_error = np.linalg.norm(FOM_model - RBF_model)/np.linalg.norm(FOM_model)
    print(f'Interpolation approximation error: {approximation_error*100} %')

def LaunchSimulations(mu_test):
    """
    This method should be parallel capable
    """
    if not os.path.exists("DataBase"):
        os.mkdir("DataBase")
    with open('ProjectParameters.json','r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    for mu in mu_test:
        model = KratosMultiphysics.Model()
        parameters_copy = UpdateProjectParameters(parameters.Clone(), mu)
        analysis_stage_class = _GetAnalysisStageClass(parameters_copy)
        simulation = CustomizeSimulation(analysis_stage_class,model,parameters_copy)
        simulation.Run()
        for process in simulation._GetListOfOutputProcesses():
            if isinstance(process, CalculateRomBasisOutputProcess):
                BasisOutputProcess = process
        np.save(f'DataBase/{mu[0]}, {mu[1]}',BasisOutputProcess._GetSnapshotsMatrix()) 

def CustomizeSimulation(cls, global_model, parameters):
    class CustomSimulation(cls):
        def __init__(self, model,project_parameters, custom_param = None):
            super().__init__(model,project_parameters)
            self.custom_param  = custom_param
        def Initialize(self):
            super().Initialize()
        def InitializeSolutionStep(self):
            super().InitializeSolutionStep()
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

def space_values(mu):

    file = f'{mu[0]}, {mu[1]}.npy'
    data_set = np.load(f'DataBase/{file}')

    if data_set.ndim != 1:
        if data_set.shape[1] > data_set.shape[0]: data_set = data_set.T
    print(data_set.shape)
    
    if data_set.ndim == 1:
        steps = 1
    else:
        steps = data_set.shape[1]

    for i in range(steps):

        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("MainModelPart")

        model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 2
        model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.0)
        model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 1.0)
        model_part.ProcessInfo.SetValue(KratosMultiphysics.STEP, 1)

        model_part.AddNodalSolutionStepVariable(CPFApp.VELOCITY_POTENTIAL)
        model_part.AddNodalSolutionStepVariable(CPFApp.AUXILIARY_VELOCITY_POTENTIAL)

        # Read the med file
        KratosMED.MedModelPartIO("model_mesh_0.med", KratosMED.MedModelPartIO.READ).ReadModelPart(model_part)# apply the elements and conditions
        params = KratosMultiphysics.Parameters("""{
            "elements_list"   : [{
                "model_part_name" : "MainModelPart.Parts_Parts_Auto1",
                "element_name"    : "Element2D3N"
            }],
            "conditions_list" : [{
                "model_part_name" : "MainModelPart.PotentialWallCondition2D_Far_field_Auto1",
                "condition_name"  : "WallCondition2D2N"
            },{
                "model_part_name" : "MainModelPart.Body2D_Body",
                "condition_name"  : "WallCondition2D2N"
            }]
        }""")

        modeler = KratosMultiphysics.CreateEntitiesFromGeometriesModeler(current_model, params)
        modeler.SetupModelPart()
        # Assign a fresh properties container to the model
        properties = model_part.CreateNewProperties(1)
        for cond in model_part.Conditions:
            cond.Properties = properties

        for elem in model_part.Elements:
            elem.Properties = properties

        space = []; values = []
        for node in model_part.Nodes:
            offset = np.where(np.arange(1,model_part.NumberOfNodes()+1, dtype=int) == node.Id)[0][0]*2
            space.append([node.X, node.Y])
            space.append([node.X, node.Y])
            if data_set.ndim == 1:
                node.SetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL, data_set[offset])
                node.SetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL, data_set[offset+1])
                values.append([data_set[offset]])
                values.append([data_set[offset+1]])
            else:
                node.SetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL, data_set[offset,i])
                node.SetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL, data_set[offset+1,i])
                values.append([data_set[offset,i]])
                values.append([data_set[offset+1,i]])

        space  = np.array(space)
        values = np.array(values)
    return space, values

if __name__ == '__main__':

    mu_train, mu_test = load_mu_parameters()

    # NumberOfMuTest = 1
    # mu_test  = random.sample(list(mu_test), NumberOfMuTest)

    #### DATA
    #################################################################################################
    parameters = np.array(mu_train)
    print(f'Parameters shape: {parameters.shape}')

    #### RBF TRAINING
    #################################################################################################
    pod = POD()
    rbf = RBF()

    db = Database()
    for mu in parameters:
        space, values = space_values(mu)
        snap = Snapshot(values=values, space=space)
        db.add(Parameter(mu), snap)

    interp = ANN([10, 10], torch.nn.Softplus(), 1000, frequency_print=200, lr=0.03)
    shift  = ANN([], torch.nn.LeakyReLU(), [1e-3, 4000], optimizer=torch.optim.Adam, frequency_print=200, l2_regularization=0,  lr=0.0023)
    nnspod = AutomaticShiftSnapshots(shift, interp, Linear(fill_value=0.0), parameter_index=2, reference_index=2, barycenter_loss=20.)
    rom    = ROM(db, pod, rbf, plugins=[nnspod])

    for _ in range(10):
        rom.fit()    # Calculate reduced space
        if rom.plugins[0].shift_network.loss_trend[-1] < 1e-3:
            break

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

    error_estimation(mu_train, mu_test)
