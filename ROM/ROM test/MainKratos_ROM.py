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


def update_project_parameters(parameters, a_case,results_name):
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(a_case[0])
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(a_case[1])
    #storing results into a results folder
    #parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString(f"Results{results_name}[{str(np.round(a_case[0]*180/math.pi,2))}] - [{str(np.round(a_case[1],2))}]")
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
        parameters = update_project_parameters(parameters, mu[i],"FOM/")
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
        parameters = update_project_parameters(parameters, mu[i],"ROM/")
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
    u,_,_,_ = RandomizedSingularValueDecomposition(COMPUTE_V=False).Calculate(RedidualsSnapshotsMatrix, 1e-12)
    simulation.GetHROM_utility().hyper_reduction_element_selector.SetUp(u)
    simulation.GetHROM_utility().hyper_reduction_element_selector.Run()
    simulation.GetHROM_utility().AppendHRomWeightsToRomParameters()
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
        parameters = update_project_parameters(parameters, mu[i],"HROM/")
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


def primal_FOM(mach):
    with open("ProjectParametersPrimalROM.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach)
    parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString('Results/FOM'+str(mach))
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

def primal_ROM(mach):
    with open("ProjectParametersPrimalROM.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach)
    parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString('Results/ROM'+str(mach))
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

def primal_HROM(mach):
    with open("ProjectParametersPrimalROM.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach)
    parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString('Results/HROM'+str(mach))
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
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('ResultsFOM')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('ResultsROM')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('ResultsHROM')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('ResultsTrain_HROM')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('ResultsHROM')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('PrimalRomParameters.json')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('ROM test.post.lst')

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

    print(":::::::::::::::::: CASE 1 :::::::::::::::::::::")
    fom_snapshots,tmfom = primal_FOM(0.3)
    setting_flags_rom_parameters(simulation_to_run = 'ROM', parameters_file_name = './PrimalRomParameters.json')
    rom_snapshots,tmrom = primal_ROM(0.3)
    setting_flags_rom_parameters(simulation_to_run = 'runHROM', parameters_file_name = './PrimalRomParameters.json')
    hrom_snapshots,tmhrom = primal_HROM(0.3)
    setting_flags_rom_parameters(simulation_to_run = 'ROM', parameters_file_name = './PrimalRomParameters.json')
    print("==========================> approximation error primal   FOM vs ROM: ",np.linalg.norm(fom_snapshots - rom_snapshots)/np.linalg.norm(fom_snapshots)*100,"%")
    print("==========================> approximation error primal  ROM vs HROM: ",np.linalg.norm(rom_snapshots - hrom_snapshots)/np.linalg.norm(rom_snapshots)*100,"%")
    print("time  FOM:", tmfom)
    print("time  ROM:", tmrom)
    print("time HROM:", tmhrom)
    print(":::::::::::::::::: CASE 2 :::::::::::::::::::::")
    fom_snapshots,tmfom = primal_FOM(0.5)
    setting_flags_rom_parameters(simulation_to_run = 'ROM', parameters_file_name = './PrimalRomParameters.json')
    rom_snapshots,tmrom = primal_ROM(0.5)
    setting_flags_rom_parameters(simulation_to_run = 'runHROM', parameters_file_name = './PrimalRomParameters.json')
    hrom_snapshots,tmhrom = primal_HROM(0.5)
    setting_flags_rom_parameters(simulation_to_run = 'ROM', parameters_file_name = './PrimalRomParameters.json')
    print("==========================> approximation error primal   FOM vs ROM: ",np.linalg.norm(fom_snapshots - rom_snapshots)/np.linalg.norm(fom_snapshots)*100,"%")
    print("==========================> approximation error primal  ROM vs HROM: ",np.linalg.norm(rom_snapshots - hrom_snapshots)/np.linalg.norm(rom_snapshots)*100,"%")
    print("time  FOM:", tmfom)
    print("time  ROM:", tmrom)
    print("time HROM:", tmhrom)
    print(":::::::::::::::::: CASE 3 :::::::::::::::::::::")
    fom_snapshots,tmfom = primal_FOM(0.7)
    setting_flags_rom_parameters(simulation_to_run = 'ROM', parameters_file_name = './PrimalRomParameters.json')
    rom_snapshots,tmrom = primal_ROM(0.7)
    setting_flags_rom_parameters(simulation_to_run = 'runHROM', parameters_file_name = './PrimalRomParameters.json')
    hrom_snapshots,tmhrom = primal_HROM(0.7)
    setting_flags_rom_parameters(simulation_to_run = 'ROM', parameters_file_name = './PrimalRomParameters.json')
    print("==========================> approximation error primal   FOM vs ROM: ",np.linalg.norm(fom_snapshots - rom_snapshots)/np.linalg.norm(fom_snapshots)*100,"%")
    print("==========================> approximation error primal  ROM vs HROM: ",np.linalg.norm(rom_snapshots - hrom_snapshots)/np.linalg.norm(rom_snapshots)*100,"%")
    print("time  FOM:", tmfom)
    print("time  ROM:", tmrom)
    print("time HROM:", tmhrom)