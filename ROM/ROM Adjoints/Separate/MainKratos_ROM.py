import KratosMultiphysics
from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis import PotentialFlowAnalysis
from KratosMultiphysics.RomApplication.rom_testing_utilities import SetUpSimulationInstance
import KratosMultiphysics.kratos_utilities
import numpy as np
from scipy.stats import qmc
import math
import json
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition



def get_multiple_params():
    sampler = qmc.Halton(d=1)
    #Angle of attack
    number_of_angles = 2
    l_angle = [-1.5 * math.pi / 180.0]
    u_angle = [ 1.5 * math.pi / 180.0]
    angle_values = qmc.scale(sampler.random(number_of_angles), l_angle, u_angle)
    #Mach infinit
    number_of_mach = 2
    l_mach = [0.15]
    u_mach = [0.25]
    mach_values = qmc.scale(sampler.random(number_of_mach), l_mach, u_mach)
    mu = []
    for i in range(number_of_angles):
        for j in range(number_of_mach):
                mu.append([angle_values[i],mach_values[j]])
    return mu



def update_project_parameters(results_name, parameters, a_case):
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(a_case[0])
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(a_case[1])
    #storing results into a results folder
    parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString('Results'+ results_name+ ' - ['+ str(np.round(a_case[0]*180/math.pi,2))+ '] - ['+ str(np.round(a_case[1],2)) +"]")
    return parameters



def multiple_params_Train_Primal_ROM(mu):
    for i in range(len(mu)):
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        print(":::::::::::::Primal solution N° ",i,":::::::::::::::")
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        with open("ProjectParametersPrimalROM.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        parameters = update_project_parameters("_Primal/" + str(i), parameters, mu[i])
        parameters["processes"]["list_other_processes"][0]["Parameters"]["file_settings"]["file_name"].SetString("Primal_h5_files/MainModelPart"+ str(i))
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

def multiple_params_Primal_ROM(mu):
    for i in range(len(mu)):
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        print("::::::::::::Primal ROM solution N° ",i,"::::::::::::")
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        with open("ProjectParametersPrimalROM.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        parameters["processes"]["list_other_processes"][0]["Parameters"]["file_settings"]["file_name"].SetString("Primal_ROM_h5_files/MainModelPart"+ str(i))
        parameters = update_project_parameters("_Primal_ROM/" + "N° "+ str(i), parameters, mu[i])
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

def multiple_params_Train_Primal_HROM(mu):
    for i in range(len(mu)):
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        print("::::::::::::Primal ROM solution N° ",i,"::::::::::::")
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        with open("ProjectParametersPrimalROM.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        parameters["processes"]["list_other_processes"][0]["Parameters"]["file_settings"]["file_name"].SetString("Primal_ROM_HROM_h5_files/MainModelPart"+ str(i))
        parameters = update_project_parameters("_Primal_ROM_HROM/" + "N° "+ str(i), parameters, mu[i])
        model = KratosMultiphysics.Model()
        simulation = SetUpSimulationInstance(model,parameters)
        simulation.Run()
        if i==0:
            RedidualsSnapshotsMatrix = simulation.GetHROM_utility()._GetResidualsProjectedMatrix() 
        else:
            RedidualsSnapshotsMatrix = np.c_[RedidualsSnapshotsMatrix, simulation.GetHROM_utility()._GetResidualsProjectedMatrix()]
    u,_,_,_ = RandomizedSingularValueDecomposition(COMPUTE_V=False).Calculate(RedidualsSnapshotsMatrix, 1e-12)
    simulation.GetHROM_utility().hyper_reduction_element_selector.SetUp(u)
    simulation.GetHROM_utility().hyper_reduction_element_selector.Run()
    simulation.GetHROM_utility().AppendHRomWeightsToRomParameters()

def multiple_params_Primal_HROM(mu):
    for i in range(len(mu)):
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        print("::::::::::::Primal HROM solution N° ",i,":::::::::::")
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        with open("ProjectParametersPrimalROM.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        parameters["processes"]["list_other_processes"][0]["Parameters"]["file_settings"]["file_name"].SetString("Primal_HROM_h5_files/MainModelPart"+ str(i))
        parameters = update_project_parameters("_Primal_HROM/" + "N° "+ str(i), parameters, mu[i])
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



def multiple_params_Train_Adjoint_ROM(mu):
    for i in range(len(mu)):
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        print("::::::::::::Adjoint solution N° ",i,":::::::::::::::")
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        with open("ProjectParametersAdjointROM.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        parameters["processes"]["list_other_processes"][0]["Parameters"]["file_settings"]["file_name"].SetString("Primal_h5_files/MainModelPart"+ str(i))
        parameters = update_project_parameters("_Adjoint/" + "N° "+ str(i), parameters, mu[i])
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

def multiple_params_Adjoint_ROM(mu):
    for i in range(len(mu)):
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        print(":::::::::::Adjoint ROM solution N° ",i,"::::::::::::")
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        with open("ProjectParametersAdjointROM.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        parameters["processes"]["list_other_processes"][0]["Parameters"]["file_settings"]["file_name"].SetString("Primal_h5_files/MainModelPart"+ str(i))
        parameters = update_project_parameters("_Adjoint_ROM/" + "N° "+ str(i), parameters, mu[i])
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

def multiple_params_Train_Adjoint_HROM(mu):
    for i in range(len(mu)):
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        print(":::::::::::Adjoint ROM solution N° ",i,"::::::::::::")
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        with open("ProjectParametersAdjointROM.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        parameters["processes"]["list_other_processes"][0]["Parameters"]["file_settings"]["file_name"].SetString("Primal_h5_files/MainModelPart"+ str(i))
        parameters = update_project_parameters("_Adjoint_ROM_HROM/" + "N° "+ str(i), parameters, mu[i])
        model = KratosMultiphysics.Model()
        simulation = SetUpSimulationInstance(model,parameters)
        simulation.Run()
        if i==0:
            RedidualsSnapshotsMatrix = simulation.GetHROM_utility()._GetResidualsProjectedMatrix() 
        else:
            RedidualsSnapshotsMatrix = np.c_[RedidualsSnapshotsMatrix, simulation.GetHROM_utility()._GetResidualsProjectedMatrix()]
    u,_,_,_ = RandomizedSingularValueDecomposition(COMPUTE_V=False).Calculate(RedidualsSnapshotsMatrix, 1e-12)
    simulation.GetHROM_utility().hyper_reduction_element_selector.SetUp(u)
    simulation.GetHROM_utility().hyper_reduction_element_selector.Run()
    simulation.GetHROM_utility().AppendHRomWeightsToRomParameters()

def multiple_params_Adjoint_HROM(mu):
    for i in range(len(mu)):
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        print(":::::::::::Adjoint HROM solution N° ",i,":::::::::::")
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        with open("ProjectParametersAdjointROM.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        parameters["processes"]["list_other_processes"][0]["Parameters"]["file_settings"]["file_name"].SetString("Primal_h5_files/MainModelPart"+ str(i))
        parameters = update_project_parameters("_Adjoint_HROM/" + "N° "+ str(i), parameters, mu[i])
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





if __name__ == "__main__":
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Results_Primal')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Results_Primal_ROM')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Results_Primal_ROM_HROM')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Results_Primal_HROM')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Results_Adjoint')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Results_Adjoint_ROM')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Results_Adjoint_ROM_HROM')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Results_Adjoint_HROM')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Primal_h5_files')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Primal_ROM_h5_files')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Primal_ROM_HROM_h5_files')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Primal_HROM_h5_files')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('PrimalRomParameters.json')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('AdjointRomParameters.json')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('ROM.post.lst')
    mu = get_multiple_params()
    primal_fom_snapshots = multiple_params_Train_Primal_ROM(mu)
    primal_rom_snapshots = multiple_params_Primal_ROM(mu)
    print("====================================> approximation error primal FOM vs ROM: ",np.linalg.norm(primal_fom_snapshots - primal_rom_snapshots)/np.linalg.norm(primal_fom_snapshots)*100,"%")
    setting_flags_rom_parameters(simulation_to_run = 'trainHROM', parameters_file_name = './PrimalRomParameters.json')
    multiple_params_Train_Primal_HROM(mu)
    setting_flags_rom_parameters(simulation_to_run = 'runHROM', parameters_file_name = './PrimalRomParameters.json')
    primal_hrom_snapshots = multiple_params_Primal_HROM(mu)
    print("====================================> approximation error primal FOM vs ROM: ",np.linalg.norm(primal_fom_snapshots - primal_rom_snapshots)/np.linalg.norm(primal_fom_snapshots)*100,"%")
    print("====================================> approximation error primal HROM vs ROM: ",np.linalg.norm(primal_rom_snapshots - primal_hrom_snapshots)/np.linalg.norm(primal_rom_snapshots)*100,"%")

    adjoint_fom_snapshots = multiple_params_Train_Adjoint_ROM(mu)
    adjoint_rom_snapshots = multiple_params_Adjoint_ROM(mu)
    print("====================================> approximation error Adjoint FOM vs ROM: ",np.linalg.norm(adjoint_fom_snapshots - adjoint_rom_snapshots)/np.linalg.norm(adjoint_fom_snapshots)*100,"%")
    setting_flags_rom_parameters(simulation_to_run = 'trainHROM', parameters_file_name = './AdjointRomParameters.json')
    multiple_params_Train_Adjoint_HROM(mu)
    setting_flags_rom_parameters(simulation_to_run = 'runHROM', parameters_file_name = './AdjointRomParameters.json')
    adjoint_hrom_snapshots = multiple_params_Adjoint_HROM(mu)
    print("====================================> approximation error primal   FOM vs ROM: ",np.linalg.norm(primal_fom_snapshots - primal_rom_snapshots)/np.linalg.norm(primal_fom_snapshots)*100,"%")
    print("====================================> approximation error primal  HROM vs ROM: ",np.linalg.norm(primal_rom_snapshots - primal_hrom_snapshots)/np.linalg.norm(primal_rom_snapshots)*100,"%")
    print("====================================> approximation error Adjoint  FOM vs ROM: ",np.linalg.norm(adjoint_fom_snapshots - adjoint_rom_snapshots)/np.linalg.norm(adjoint_fom_snapshots)*100,"%")
    print("====================================> approximation error Adjoint HROM vs ROM: ",np.linalg.norm(adjoint_rom_snapshots - adjoint_hrom_snapshots)/np.linalg.norm(adjoint_rom_snapshots)*100,"%")