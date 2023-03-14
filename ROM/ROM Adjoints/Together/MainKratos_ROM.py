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



def multiple_params_TrainROM(mu):
    for i in range(len(mu)):
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        print(":::::::::::::Primal solution N° ",i,":::::::::::::::")
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        with open("ProjectParametersPrimalROM.json",'r') as primal_parameter_file:
            primal_parameters = KratosMultiphysics.Parameters(primal_parameter_file.read())
        primal_parameters = update_project_parameters("_Primal/" + str(i), primal_parameters, mu[i])
        primal_parameters["processes"]["list_other_processes"][0]["Parameters"]["file_settings"]["file_name"].SetString("Primal_h5_files/MainModelPart"+ str(i))
        primal_model = KratosMultiphysics.Model()
        primal_simulation = PotentialFlowAnalysis(primal_model,primal_parameters)
        primal_simulation.Run()
        for primal_process in primal_simulation._GetListOfOutputProcesses():
            if type(primal_process) == KratosMultiphysics.RomApplication.calculate_rom_basis_output_process.CalculateRomBasisOutputProcess:
                Primal_BasisOutputProcess = primal_process
        if i==0:
            Primal_SnapshotsMatrix = Primal_BasisOutputProcess._GetSnapshotsMatrix()
        else:
            Primal_SnapshotsMatrix = np.c_[Primal_SnapshotsMatrix, Primal_BasisOutputProcess._GetSnapshotsMatrix()]
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        print("::::::::::::Adjoint solution N° ",i,":::::::::::::::")
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        with open("ProjectParametersAdjointROM.json",'r') as adjoint_parameter_file:
            adjoint_parameters = KratosMultiphysics.Parameters(adjoint_parameter_file.read())
        adjoint_parameters = update_project_parameters("_Adjoints/" + "N° "+ str(i), adjoint_parameters, mu[i])
        adjoint_parameters["processes"]["list_other_processes"][0]["Parameters"]["file_settings"]["file_name"].SetString("Primal_h5_files/MainModelPart"+ str(i))
        adjoint_model = KratosMultiphysics.Model()
        adjoint_simulation = PotentialFlowAnalysis(adjoint_model,adjoint_parameters)
        adjoint_simulation.Run()
        for adjoint_process in adjoint_simulation._GetListOfOutputProcesses():
            if type(adjoint_process) == KratosMultiphysics.RomApplication.calculate_rom_basis_output_process.CalculateRomBasisOutputProcess:
                Adjoint_BasisOutputProcess = adjoint_process
        if i==0:
            Adjoint_SnapshotsMatrix = Adjoint_BasisOutputProcess._GetSnapshotsMatrix()
        else:
            Adjoint_SnapshotsMatrix = np.c_[Adjoint_SnapshotsMatrix, Adjoint_BasisOutputProcess._GetSnapshotsMatrix()]
    Primal_BasisOutputProcess._PrintRomBasis(Primal_SnapshotsMatrix)
    Adjoint_BasisOutputProcess._PrintRomBasis(Adjoint_SnapshotsMatrix)
    return (Primal_SnapshotsMatrix, Adjoint_SnapshotsMatrix)



def multiple_params_TrainHROM(mu):
    for i in range(len(mu)):
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        print("::::::::::::Primal ROM solution N° ",i,"::::::::::::")
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        with open("ProjectParametersPrimalROM.json",'r') as primal_parameter_file:
            primal_parameters = KratosMultiphysics.Parameters(primal_parameter_file.read())
        primal_parameters = update_project_parameters("_ROM_Primal/" + "N° "+ str(i), primal_parameters, mu[i])
        primal_parameters["processes"]["list_other_processes"][0]["Parameters"]["file_settings"]["file_name"].SetString("Primal_ROM_h5_files/MainModelPart"+ str(i))
        primal_model = KratosMultiphysics.Model()
        primal_simulation = SetUpSimulationInstance(primal_model,primal_parameters)
        primal_simulation.Run()
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        print(":::::::::::Adjoint ROM solution N° ",i,"::::::::::::")
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        with open("ProjectParametersAdjointROM.json",'r') as Adjoint_parameter_file:
            Adjoint_parameters = KratosMultiphysics.Parameters(Adjoint_parameter_file.read())
        Adjoint_parameters = update_project_parameters("_ROM_Adjoint/" + "N° "+ str(i), Adjoint_parameters, mu[i])
        Adjoint_parameters["processes"]["list_other_processes"][0]["Parameters"]["file_settings"]["file_name"].SetString("Primal_h5_files/MainModelPart"+ str(i))
        Adjoint_model = KratosMultiphysics.Model()
        Adjoint_simulation = SetUpSimulationInstance(Adjoint_model,Adjoint_parameters)
        Adjoint_simulation.Run()
        for primal_process in primal_simulation._GetListOfOutputProcesses():
            if type(primal_process) == KratosMultiphysics.RomApplication.calculate_rom_basis_output_process.CalculateRomBasisOutputProcess:
                primal_BasisOutputProcess = primal_process
        for adjoint_process in Adjoint_simulation._GetListOfOutputProcesses():
            if type(adjoint_process) == KratosMultiphysics.RomApplication.calculate_rom_basis_output_process.CalculateRomBasisOutputProcess:
                Adjoint_BasisOutputProcess = adjoint_process
        if i==0:
            Primal_SnapshotsMatrix = primal_BasisOutputProcess._GetSnapshotsMatrix()
            primal_RedidualsSnapshotsMatrix = primal_simulation.GetHROM_utility()._GetResidualsProjectedMatrix()
            Adjoint_SnapshotsMatrix = Adjoint_BasisOutputProcess._GetSnapshotsMatrix()
            Adjoint_RedidualsSnapshotsMatrix = Adjoint_simulation.GetHROM_utility()._GetResidualsProjectedMatrix()
        else:
            Primal_SnapshotsMatrix = np.c_[Primal_SnapshotsMatrix, primal_BasisOutputProcess._GetSnapshotsMatrix()]
            primal_RedidualsSnapshotsMatrix = np.c_[primal_RedidualsSnapshotsMatrix, primal_simulation.GetHROM_utility()._GetResidualsProjectedMatrix()]
            Adjoint_SnapshotsMatrix = np.c_[Adjoint_SnapshotsMatrix, Adjoint_BasisOutputProcess._GetSnapshotsMatrix()]
            Adjoint_RedidualsSnapshotsMatrix = np.c_[Adjoint_RedidualsSnapshotsMatrix, Adjoint_simulation.GetHROM_utility()._GetResidualsProjectedMatrix()]
    u,_,_,_ = RandomizedSingularValueDecomposition(COMPUTE_V=False).Calculate(primal_RedidualsSnapshotsMatrix, 1e-12)
    u,_,_,_ = RandomizedSingularValueDecomposition(COMPUTE_V=False).Calculate(Adjoint_RedidualsSnapshotsMatrix, 1e-12)
    primal_simulation.GetHROM_utility().hyper_reduction_element_selector.SetUp(u)
    primal_simulation.GetHROM_utility().hyper_reduction_element_selector.Run()
    primal_simulation.GetHROM_utility().AppendHRomWeightsToRomParameters()
    Adjoint_simulation.GetHROM_utility().hyper_reduction_element_selector.SetUp(u)
    Adjoint_simulation.GetHROM_utility().hyper_reduction_element_selector.Run()
    Adjoint_simulation.GetHROM_utility().AppendHRomWeightsToRomParameters()
    return (Primal_SnapshotsMatrix, Adjoint_SnapshotsMatrix)



def multiple_params_HROM(mu):
    for i in range(len(mu)):
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        print("::::::::::::Primal HROM solution N° ",i,":::::::::::")
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        with open("ProjectParametersPrimalROM.json",'r') as primal_parameter_file:
            primal_parameters = KratosMultiphysics.Parameters(primal_parameter_file.read())
        primal_parameters = update_project_parameters("_HROM_Primal/" + "N° "+ str(i), primal_parameters, mu[i])
        primal_parameters["processes"]["list_other_processes"][0]["Parameters"]["file_settings"]["file_name"].SetString("Primal_HROM_h5_files/MainModelPart"+ str(i))
        primal_model = KratosMultiphysics.Model()
        primal_simulation = SetUpSimulationInstance(primal_model,primal_parameters)
        primal_simulation.Run()
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        print(":::::::::::Adjoint HROM solution N° ",i,":::::::::::")
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        with open("ProjectParametersAdjointROM.json",'r') as Adjoint_parameter_file:
            Adjoint_parameters = KratosMultiphysics.Parameters(Adjoint_parameter_file.read())
        Adjoint_parameters = update_project_parameters("_HROM_Adjoint/" + "N° "+ str(i), Adjoint_parameters, mu[i])
        Adjoint_parameters["processes"]["list_other_processes"][0]["Parameters"]["file_settings"]["file_name"].SetString("Primal_h5_files/MainModelPart"+ str(i))
        Adjoint_model = KratosMultiphysics.Model()
        Adjoint_simulation = SetUpSimulationInstance(Adjoint_model,Adjoint_parameters)
        Adjoint_simulation.Run()
        for process in primal_simulation._GetListOfOutputProcesses():
            if type(process) == KratosMultiphysics.RomApplication.calculate_rom_basis_output_process.CalculateRomBasisOutputProcess:
                primal_BasisOutputProcess = process
        for process in Adjoint_simulation._GetListOfOutputProcesses():
            if type(process) == KratosMultiphysics.RomApplication.calculate_rom_basis_output_process.CalculateRomBasisOutputProcess:
                Adjoint_BasisOutputProcess = process
        if i==0:
            Primal_SnapshotsMatrix = primal_BasisOutputProcess._GetSnapshotsMatrix()
            Adjoint_SnapshotsMatrix = Adjoint_BasisOutputProcess._GetSnapshotsMatrix()
        else:
            Primal_SnapshotsMatrix = np.c_[Primal_SnapshotsMatrix, primal_BasisOutputProcess._GetSnapshotsMatrix()]
            Adjoint_SnapshotsMatrix = np.c_[Adjoint_SnapshotsMatrix, Adjoint_BasisOutputProcess._GetSnapshotsMatrix()]
    return (Primal_SnapshotsMatrix, Adjoint_SnapshotsMatrix)



def setting_flags_rom_parameters(simulation_to_run = 'ROM'):
    parameters_file_name = './PrimalRomParameters.json'
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
    parameters_file_name = './AdjointRomParameters.json'
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
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Results_Adjoints')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Results_ROM_Primal')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Results_ROM_Adjoint')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Results_HROM_Primal')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Results_HROM_Adjoint')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Primal_h5_files')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Primal_ROM_h5_files')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Primal_HROM_h5_files')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('PrimalRomParameters.json')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('AdjointRomParameters.json')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('ROM2.post.lst')
    mu = get_multiple_params()
    fom_snapshots = multiple_params_TrainROM(mu)
    setting_flags_rom_parameters(simulation_to_run = 'trainHROM')
    rom_snapshots = multiple_params_TrainHROM(mu)
    setting_flags_rom_parameters(simulation_to_run = 'runHROM')
    hrom_snapshots = multiple_params_HROM(mu)
    setting_flags_rom_parameters(simulation_to_run = 'ROM')
    print(":::::::::::::::::::::::: Primal cases ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print("==============================> approximation error FOM vs ROM: ",np.linalg.norm(fom_snapshots[0] - rom_snapshots[0])/np.linalg.norm(fom_snapshots[0])*100,"%")
    print("==============================> approximation error HROM vs ROM: ",np.linalg.norm(rom_snapshots[0] - hrom_snapshots[0])/np.linalg.norm(rom_snapshots[0])*100,"%")
    print("::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print(":::::::::::::::::::::::: Adjoint cases :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    print("==============================> approximation error FOM vs ROM: ",np.linalg.norm(fom_snapshots[1] - rom_snapshots[1])/np.linalg.norm(fom_snapshots[1])*100,"%")
    print("==============================> approximation error HROM vs ROM: ",np.linalg.norm(rom_snapshots[1] - hrom_snapshots[1])/np.linalg.norm(rom_snapshots[1])*100,"%")

