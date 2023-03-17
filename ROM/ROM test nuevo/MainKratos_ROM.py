import KratosMultiphysics
from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis import PotentialFlowAnalysis
from KratosMultiphysics.RomApplication.rom_testing_utilities import SetUpSimulationInstance
import KratosMultiphysics.kratos_utilities
import numpy as np
from scipy.stats import qmc
import matplotlib.pyplot as plt
import math
import time
import json
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition
from KratosMultiphysics.RomApplication.empirical_cubature_method import EmpiricalCubatureMethod
import KratosMultiphysics.RomApplication as KratosROM

def Parallel_ECM_Recursive(arr,block_len,global_ids,global_weights,final=False):
    # Arr is a dislib array with the chunks coinciding with the row partitions we want

    number_of_rows=0
    for block in arr:
        number_of_rows+=block.shape[0]

    if final:
        print('STARTED FINAL ECM')
        print('Chunks in dslibarray:', len(arr))
        ecm_mode='final'
        if len(arr)>1:
            print('Non single-chunk, rechunking.')
            arr=[np.concatenate(arr, axis=0)]
            # arr=arr.rechunk(arr.shape)
        block_len=number_of_rows
    else:
        print('STARTED INTERMEDIATE ECM')
        ecm_mode='intermediate'

    # Creating lists to store ids and weights
    ids_global_store = []
    ids_local_store= []
    weights_store = []
    
    for j in range(len(arr)):
        # We get the part of the global ids and weights corresponding to the elements in the chunk.
        chunk_global_ids = global_ids[list(range(block_len*j,min(block_len*(j+1),number_of_rows)))]
        chunk_global_weights = global_weights[list(range(block_len*j,min(block_len*(j+1),number_of_rows)))]
        # We run the ECM algorithm and get the list of chosen elements in both local and global notation, and their global weights
        ids,l_ids,w = GetElementsOfPartition(arr[j],chunk_global_ids,chunk_global_weights,block_len,j,ecm_mode)
        ids_global_store.append(ids)
        ids_local_store.append(l_ids)
        weights_store.append(w)

    # Synchronize the ids and weights lists
    for i in range(len(ids_global_store)):
        if i==0:
            temp_global_ids = ids_global_store[i]
            temp_local_ids = ids_local_store[i]
            temp_global_weights = weights_store[i]
        else:
            temp_global_ids = np.r_[temp_global_ids,ids_global_store[i]]
            temp_local_ids = np.r_[temp_local_ids,ids_local_store[i]]
            temp_global_weights = np.r_[temp_global_weights,weights_store[i]]

    global_ids = temp_global_ids
    local_ids = temp_local_ids
    global_weights = temp_global_weights

    # We return the rows of the array corresponding to the chosen elements. Also the global ids and weights.
    # This dislib array will still have the the same chunk size, just with less chunks
    result_mat=np.concatenate(arr, axis=0)
    result_mat=result_mat[local_ids]
    arr=[]
    k=0
    while k*block_len < len(result_mat):
        arr.append(result_mat[k*block_len:min(len(result_mat)-1,(k+1)*block_len)])
        k+=1
    return arr, global_ids, global_weights


def GetElementsOfPartition(np_array, global_ids, global_weights, block_len, block_num, title):

    projected_residuals_matrix = np_array * global_weights[:, np.newaxis] #try making the weights and indexes ds arrays?
    
    #u,_,_ = truncated_svd(projected_residuals_matrix, 0)  #numpy "exact" version
    if title == 'intermediate':
        #u, _, _, _ = RandomizedSingularValueDecomposition().Calculate(projected_residuals_matrix) #randomized version with machine precision
        M, N = projected_residuals_matrix.shape 
        u,s,_ = np.linalg.svd(projected_residuals_matrix)
        k = get_number_of_singular_values_for_given_tolerance(M, N, s, 0)
        u = u[:,:k]
        constrain_sum_of_weights = False # setting it to "True" worsens the approximation. Need to implement the orthogonal complement rather and not the row of 1's is implemented
    else:
        #u, _, _, _ = RandomizedSingularValueDecomposition().Calculate(projected_residuals_matrix,1e-6) #randomized version with user-defined tolerance
        M, N = projected_residuals_matrix.shape 
        u,s,_ = np.linalg.svd(projected_residuals_matrix)
        k = get_number_of_singular_values_for_given_tolerance(M, N, s, 1e-12)
        u = u[:,:k]        
        constrain_sum_of_weights = False # setting it to "True" worsens the approximation. Need to implement the orthogonal complement rather and not the row of 1's is implemented

    ElementSelector = EmpiricalCubatureMethod()
    ElementSelector.SetUp( u, constrain_sum_of_weights)
    ElementSelector.Initialize()
    ElementSelector.Calculate()
    local_ids = np.squeeze(ElementSelector.z)
    weights = np.squeeze(ElementSelector.w)

    """ indexes_2 = np.argsort(local_ids)
    return local_ids[indexes_2], weights[indexes_2]) """
    indexes_2 = np.argsort(local_ids) #this is necessary, since dislib cannot return un-ordered indexes

    return global_ids[local_ids[indexes_2]], local_ids[indexes_2]+block_len*block_num, weights[indexes_2]*global_weights[local_ids[indexes_2]]

def Initialize_ECM_Lists(arr): # Generate first lists of ids and weights
    # number_of_rows = arr.shape[0]
    number_of_rows=0
    for block in arr:
        number_of_rows+=block.shape[0]
    global_ids = np.array(range(number_of_rows))
    global_weights = np.ones(len(global_ids))

    return global_ids, global_weights

def AppendHRomWeightsToRomParameters(hrom_training_utility, z, w):
        n_elements = hrom_training_utility.solver.GetComputingModelPart().NumberOfElements()

        # Create dictionary with HROM weights
        hrom_weights = {}
        hrom_weights["Elements"] = {}
        hrom_weights["Conditions"] = {}

        if type(z)==np.int64 or type(z)==np.int32:
            # Only one element found !
            if z <= n_elements-1:
                hrom_weights["Elements"][int(z)] = float(w)
            else:
                hrom_weights["Conditions"][int(z)-n_elements] = float(w)
        else:
            # Many elements found
            for j in range (0,len(z)):
                if z[j] <=  n_elements -1:
                    hrom_weights["Elements"][int(z[j])] = float(w[j])
                else:
                    hrom_weights["Conditions"][int(z[j])-n_elements] = float(w[j])

        #TODO: Make this optional
        # If required, keep at least one condition per submodelpart
        # This might be required by those BCs involving the faces (e.g. slip BCs)
        include_minimum_condition = False
        if include_minimum_condition:
            # Get the HROM conditions to be added
            minimum_conditions = KratosROM.RomAuxiliaryUtilities.GetHRomMinimumConditionsIds(
                hrom_training_utility.solver.GetComputingModelPart().GetRootModelPart(), #TODO: I think this one should be the root
                hrom_weights["Conditions"])

            # Add the selected conditions to the conditions dict with a null weight
            for cond_id in minimum_conditions:
                hrom_weights["Conditions"][cond_id] = 0.0

        #TODO: Make this optional
        # If required, add the HROM conditions parent elements
        # Note that we add these with zero weight so their future assembly will have no effect
        include_condition_parents = False
        if include_condition_parents:
            # Get the HROM condition parents from the current HROM weights
            missing_condition_parents = KratosROM.RomAuxiliaryUtilities.GetHRomConditionParentsIds(
                hrom_training_utility.solver.GetComputingModelPart().GetRootModelPart(), #TODO: I think this one should be the root
                hrom_weights)

            # Add the missing parents to the elements dict with a null weight
            for parent_id in missing_condition_parents:
                hrom_weights["Elements"][parent_id] = 0.0

        # Append weights to RomParameters.json
        # We first parse the current RomParameters.json to then append and edit the data
        with open(hrom_training_utility.rom_parameters_file_name + '.json','r') as f:
            updated_rom_parameters = json.load(f)
            #FIXME: I don't really like to automatically change things without the user realizing...
            #FIXME: However, this leaves the settings ready for the HROM postprocess... something that is cool
            #FIXME: Decide about this
            # updated_rom_parameters["train_hrom"] = False
            # updated_rom_parameters["run_hrom"] = True
            updated_rom_parameters["elements_and_weights"] = hrom_weights #TODO: Rename elements_and_weights to hrom_weights

        with open(hrom_training_utility.rom_parameters_file_name + '.json','w') as f:
            json.dump(updated_rom_parameters, f, indent = 4)

        if hrom_training_utility.echo_level > 0 : KratosMultiphysics.Logger.PrintInfo("HRomTrainingUtility","\'RomParameters.json\' file updated with HROM weights.")


def get_multiple_params():
    plot_values = True
    number_of_values = 15
    sampler = qmc.Halton(d=2)
    sample = sampler.random(number_of_values)
    #Angle of attack
    l_angle = [-0.01 * math.pi / 180.0]
    u_angle = [ 0.01 * math.pi / 180.0]
    #Mach infinit
    l_mach = [0.03]
    u_mach = [0.6]
    mu = []
    values = qmc.scale(sample, [l_angle[0],l_mach[0]], [u_angle[0],u_mach[0]])
    for i in range(number_of_values):
        mu.append([values[i,0], values[i,1]])
    if plot_values:
        for i in range(len(values)):
            plt.plot(values[i,1], np.round(values[i,0]*180/math.pi,1)+5, 'bs')
        plt.ylabel('Alpha')
        plt.xlabel('Mach')
        plt.show()
    return mu

def get_multiple_params_angle():
    plot_values = True
    number_of_values = 15
    sampler = qmc.Halton(d=1)
    sample = sampler.random(number_of_values)
    #Angle of attack
    l_angle = [-6.0 * math.pi / 180.0]
    u_angle = [ 1.0 * math.pi / 180.0]
    mu = []
    values = qmc.scale(sample, [l_angle[0]], [u_angle[0]])
    for i in range(number_of_values):
        mu.append([values[i], 0.3])
    if plot_values:
        for i in range(len(values)):
            plt.plot(np.round(values[i]*180/math.pi,1)+5, 0.3, 'bs')
        plt.xlabel('Alpha')
        plt.ylabel('Mach')
        plt.show()
    return mu

def get_multiple_params_mach():
    plot_values = False
    number_of_values = 4
    sampler = qmc.Halton(d=1)
    sample = sampler.random(number_of_values)
    #Mach infinit
    l_mach = [0.03]
    u_mach = [0.6]
    mu = []
    values = qmc.scale(sample, [l_mach[0]], [u_mach[0]])
    for i in range(number_of_values):
        mu.append([0.0, values[i]])
    if plot_values:
        for i in range(len(values)):
            plt.plot(values[i], 5.0, 'bs')
        plt.ylabel('Alpha')
        plt.xlabel('Mach')
        plt.show()
    return mu

def get_number_of_singular_values_for_given_tolerance(M, N, s, epsilon):
    dimMATRIX = max(M,N)
    tol = dimMATRIX*np.finfo(float).eps*max(s)/2
    R = np.sum(s > tol)  # Definition of numerical rank
    if epsilon == 0:
        K = R
    else:
        SingVsq = np.multiply(s,s)
        SingVsq.sort()
        normEf2 = np.sqrt(np.cumsum(SingVsq))
        epsilon = epsilon*normEf2[-1] #relative tolerance
        T = (sum(normEf2<epsilon))
        K = len(s)-T
    K = min(R,K)
    return K


def update_project_parameters(parameters, a_case):
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(a_case[0])
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(a_case[1])
    return parameters


def multiple_params_Train_Primal_ROM(mu):
    with open("ProjectParametersPrimalROM.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    for i in range(len(mu)):
        if parameters["solver_settings"].Has("element_replace_settings"):
            parameters["solver_settings"].RemoveValue("element_replace_settings")
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        print(":::::::::::::Primal solution N° ",i,":::::::::::::::")
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        parameters = update_project_parameters(parameters, mu[i])
        model = KratosMultiphysics.Model()
        simulation = PotentialFlowAnalysis(model, parameters)
        simulation.Run()
        for process in simulation._GetListOfOutputProcesses():
            if type(process) == KratosMultiphysics.RomApplication.calculate_rom_basis_output_process.CalculateRomBasisOutputProcess:
                BasisOutputProcess = process
        if i==0:
            SnapshotsMatrix = BasisOutputProcess._GetSnapshotsMatrix()
        else:
            SnapshotsMatrix = np.c_[SnapshotsMatrix, BasisOutputProcess._GetSnapshotsMatrix()]
    BasisOutputProcess._PrintRomBasis(SnapshotsMatrix)
    return SnapshotsMatrix

def multiple_params_Train_Primal_HROM(mu):
    with open("ProjectParametersPrimalROM.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    for i in range(len(mu)):
        if parameters["solver_settings"].Has("element_replace_settings"):
            parameters["solver_settings"].RemoveValue("element_replace_settings")
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        print("::::::::::::Train HROM solution N° ",i,"::::::::::::")
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        with open("ProjectParametersPrimalROM.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        parameters = update_project_parameters(parameters, mu[i])
        model = KratosMultiphysics.Model()
        simulation = SetUpSimulationInstance(model,parameters)
        simulation.Run()
        for process in simulation._GetListOfOutputProcesses():
            if type(process) == KratosMultiphysics.RomApplication.calculate_rom_basis_output_process.CalculateRomBasisOutputProcess:
                BasisOutputProcess = process
        if i==0:
            SnapshotsMatrix = BasisOutputProcess._GetSnapshotsMatrix()
            RedidualsSnapshotsMatrix = simulation.GetHROM_utility()._GetResidualsProjectedMatrix() 
        else:
            SnapshotsMatrix = np.c_[SnapshotsMatrix, BasisOutputProcess._GetSnapshotsMatrix()]
            RedidualsSnapshotsMatrix = np.c_[RedidualsSnapshotsMatrix, simulation.GetHROM_utility()._GetResidualsProjectedMatrix()]
    # u,_,_,_ = RandomizedSingularValueDecomposition(COMPUTE_V=False).Calculate(RedidualsSnapshotsMatrix, 1e-6)

    NumberOfPartitions = 6
    expected_shape=RedidualsSnapshotsMatrix.shape
    desired_block_size = (math.ceil(expected_shape[0]/NumberOfPartitions),expected_shape[1])

    # desired_block_size = (math.ceil(expected_shape[0]/NumberOfPartitions)-1,expected_shape[1])
    print('Total array size:', expected_shape)
    print('Applying chunk size:', desired_block_size)
    # arr = load_blocks_rechunk([RedidualsSnapshotsMatrix], shape = expected_shape, block_size = expected_shape, new_block_size = desired_block_size)
    
    block_len=desired_block_size[0]
    j=0
    arr=[]
    while j*desired_block_size[0] < len(RedidualsSnapshotsMatrix):
        arr.append(RedidualsSnapshotsMatrix[j*block_len:min(len(RedidualsSnapshotsMatrix),(j+1)*block_len)])
        j+=1

    print('BEFORE REDUCTION')
    print('Chunks in dslibarray:', len(arr))

    ecm_iterations = 16
    z,w = Initialize_ECM_Lists(arr)
    for i in range(ecm_iterations):
        print(i)
        print(ecm_iterations)
        if i < ecm_iterations-1:
            arr,z,w = Parallel_ECM_Recursive(arr,block_len,z,w)
            print('Chunks in dslibarray:', len(arr))
            print('Global_ids shape:', len(z))
        else:
            _,z,w = Parallel_ECM_Recursive(arr,block_len,z,w,final=True)
            print('Global_ids shape:', len(z))

    AppendHRomWeightsToRomParameters(simulation.GetHROM_utility(), z, w)
    # SavingElementsAndWeights("",len(model.GetModelPart("MainModelPart").GetElements()),z,w)

    #M, N = RedidualsSnapshotsMatrix.shape
    #u,s,_ = np.linalg.svd(RedidualsSnapshotsMatrix)
    #k = get_number_of_singular_values_for_given_tolerance(M, N, s, 1e-12)
    #u = u[:,:k]

    # simulation.GetHROM_utility().hyper_reduction_element_selector.SetUp(u)
    # simulation.GetHROM_utility().hyper_reduction_element_selector.Run()
    # simulation.GetHROM_utility().AppendHRomWeightsToRomParameters()
    return SnapshotsMatrix

def multiple_params_Primal_HROM(mu):
    with open("ProjectParametersPrimalROM.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    for i in range(len(mu)):
        if parameters["solver_settings"].Has("element_replace_settings"):
            parameters["solver_settings"].RemoveValue("element_replace_settings")
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        print("::::::::::::Primal HROM solution N° ",i,":::::::::::")
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        parameters = update_project_parameters(parameters, mu[i])
        model = KratosMultiphysics.Model()
        simulation = SetUpSimulationInstance(model,parameters) 
        simulation.Run()
        for process in simulation._GetListOfOutputProcesses():
            if type(process) == KratosMultiphysics.RomApplication.calculate_rom_basis_output_process.CalculateRomBasisOutputProcess:
                BasisOutputProcess = process
        if i==0:
            SnapshotsMatrix = BasisOutputProcess._GetSnapshotsMatrix()
        else:
            SnapshotsMatrix = np.c_[SnapshotsMatrix, BasisOutputProcess._GetSnapshotsMatrix()]
    return SnapshotsMatrix

def setting_flags_rom_parameters(simulation_to_run = 'ROM', parameters_file_name = './PrimalRomParameters.json'):
    with open(parameters_file_name, 'r+') as parameter_file:
        f=json.load(parameter_file)
        if simulation_to_run=='ROM':
            f['train_hrom']=False
            f['run_hrom']=False
        elif simulation_to_run=='trainHROM':
            f['train_hrom']=True
            f['run_hrom']=False
        elif simulation_to_run=='runHROM':
            f['train_hrom']=False
            f['run_hrom']=True
        else:
            print('Unknown operation. Add new rule!')
        parameter_file.seek(0)
        json.dump(f,parameter_file,indent=4)
        parameter_file.truncate()


def primal_FOM():
    with open("ProjectParametersPrimalROM.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString('Results/FOM')
    model = KratosMultiphysics.Model()
    simulation = PotentialFlowAnalysis(model, parameters)
    start=time.time()
    simulation.Run()
    end=time.time()  
    for process in simulation._GetListOfOutputProcesses():
        if type(process) == KratosMultiphysics.RomApplication.calculate_rom_basis_output_process.CalculateRomBasisOutputProcess:
            BasisOutputProcess = process
    SnapshotsMatrix = BasisOutputProcess._GetSnapshotsMatrix()
    tm = end - start
    return SnapshotsMatrix,tm

def primal_ROM():
    with open("ProjectParametersPrimalROM.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString('Results/ROM')
    model = KratosMultiphysics.Model()
    simulation = SetUpSimulationInstance(model,parameters) 
    start=time.time()
    simulation.Run()
    end=time.time()  
    for process in simulation._GetListOfOutputProcesses():
        if type(process) == KratosMultiphysics.RomApplication.calculate_rom_basis_output_process.CalculateRomBasisOutputProcess:
            BasisOutputProcess = process
    SnapshotsMatrix = BasisOutputProcess._GetSnapshotsMatrix()
    tm = end - start
    return SnapshotsMatrix,tm

def primal_HROM():
    with open("ProjectParametersPrimalROM.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString('Results/HROM')
    model = KratosMultiphysics.Model()
    simulation = SetUpSimulationInstance(model,parameters) 
    start=time.time()
    simulation.Run()
    end=time.time()  
    for process in simulation._GetListOfOutputProcesses():
        if type(process) == KratosMultiphysics.RomApplication.calculate_rom_basis_output_process.CalculateRomBasisOutputProcess:
            BasisOutputProcess = process
    SnapshotsMatrix = BasisOutputProcess._GetSnapshotsMatrix()
    tm = end - start
    return SnapshotsMatrix,tm




if __name__ == "__main__":
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Results')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('PrimalRomParameters.json')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('ROM test nuevo.post.lst')

    #mu = get_multiple_params()
    #mu = get_multiple_params_angle()
    mu = get_multiple_params_mach()
     
    
    primal_fom_snapshots = multiple_params_Train_Primal_ROM(mu)


    setting_flags_rom_parameters(simulation_to_run = 'trainHROM', parameters_file_name = './PrimalRomParameters.json')
    primal_rom_snapshots = multiple_params_Train_Primal_HROM(mu)
    setting_flags_rom_parameters(simulation_to_run = 'runHROM', parameters_file_name = './PrimalRomParameters.json')
    primal_hrom_snapshots = multiple_params_Primal_HROM(mu)
    setting_flags_rom_parameters(simulation_to_run = 'ROM', parameters_file_name = './PrimalRomParameters.json')
    print("==========================> approximation error primal   FOM vs ROM: ",np.linalg.norm(primal_fom_snapshots - primal_rom_snapshots)/np.linalg.norm(primal_fom_snapshots)*100,"%")
    print("==========================> approximation error primal  ROM vs HROM: ",np.linalg.norm(primal_rom_snapshots - primal_hrom_snapshots)/np.linalg.norm(primal_rom_snapshots)*100,"%")

    #fom_snapshots,tmfom = primal_FOM()
    #setting_flags_rom_parameters(simulation_to_run = 'ROM', parameters_file_name = './PrimalRomParameters.json')
    #rom_snapshots,tmrom = primal_ROM()
    #setting_flags_rom_parameters(simulation_to_run = 'runHROM', parameters_file_name = './PrimalRomParameters.json')
    #hrom_snapshots,tmhrom = primal_HROM()
    #setting_flags_rom_parameters(simulation_to_run = 'ROM', parameters_file_name = './PrimalRomParameters.json')
    #print("==========================> approximation error primal   FOM vs ROM: ",np.linalg.norm(fom_snapshots - rom_snapshots)/np.linalg.norm(fom_snapshots)*100,"%")
    #print("==========================> approximation error primal  ROM vs HROM: ",np.linalg.norm(rom_snapshots - hrom_snapshots)/np.linalg.norm(rom_snapshots)*100,"%")
    #print("time  FOM:", tmfom)
    #print("time  ROM:", tmrom)
    #print("time HROM:", tmhrom)