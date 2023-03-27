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
from KratosMultiphysics.gid_output_process import GiDOutputProcess
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition


def get_multiple_params():
    # plot_values = True
    # number_of_values = 4
    # sampler = qmc.Halton(d=1)
    # sample = sampler.random(number_of_values)
    # #Angle of attack
    # l_angle = [-6]
    # u_angle = [ 1]
    # #Mach infinit
    # l_mach = [0.03]
    # u_mach = [0.6]
    # mu = []
    # params = np.zeros((15))
    # for i in range(len(params)):
    #     params[i] = l_angle[0] + i * 0.5
        
    # for i in range(len(params)): 
    #     sample = sampler.random(number_of_values)
    #     values = qmc.scale(sample, [l_mach[0]], [u_mach[0]])
    #     mu.append([params[i] * math.pi / 180.0, l_mach[0]]) 
    #     mu.append([params[i] * math.pi / 180.0, u_mach[0]]) 
    #     for j in range(number_of_values):
    #         mu.append([params[i] * math.pi / 180.0, values[j]])  
        
    # for i in range(len(params)*(number_of_values+2)):
    #     if plot_values:
    #         plt.plot(mu[i][1], np.round(mu[i][0]*180/math.pi,2) + 5, 'bs')

    # plt.ylabel('Alpha')
    # plt.xlabel('Mach')
    # plt.grid(True)
    # plt.show()

    plot_values = True
    number_of_values = 16
    sampler = qmc.Halton(d=2)
    sample = sampler.random(number_of_values)
    #Angle of attack
    l_angle = [-6.0 * math.pi / 180.0]
    u_angle = [ 1.0 * math.pi / 180.0]
    #Mach infinit
    l_mach = [0.03]
    u_mach = [0.6]
    mu = []
    values = qmc.scale(sample, [l_angle[0],l_mach[0]], [u_angle[0],u_mach[0]])
    values[0,0] = l_angle[0]
    values[0,1] = l_mach[0]
    values[1,0] = l_angle[0]
    values[1,1] = u_mach[0]
    values[number_of_values-1,0] = u_angle[0]
    values[number_of_values-1,1] = u_mach[0]
    values[number_of_values-2,0] = u_angle[0]
    values[number_of_values-2,1] = l_mach[0]
    for i in range(number_of_values):
        mu.append([values[i,0], values[i,1]])
    if plot_values:
        for i in range(len(values)):
            plt.plot(values[i,1], np.round(values[i,0]*180/math.pi,1) + 5, 'bs')
        plt.ylabel('Alpha')
        plt.xlabel('Mach')
        plt.grid(True)
        plt.show()
    return mu

def get_multiple_params_angle():
    plot_values = True
    number_of_values = 4
    sampler = qmc.Halton(d=1)
    sample = sampler.random(number_of_values)
    #Angle of attack
    l_angle = [-5.95 * math.pi / 180.0] # -1 
    u_angle = [ 0.95 * math.pi / 180.0] #  6 
    mu = []
    values = qmc.scale(sample, [l_angle[0]], [u_angle[0]])

    mu.append([-6.0 * math.pi / 180.0, 0.3])
    # mu.append([-5.0 * math.pi / 180.0, 0.3])
    # mu.append([-4.0 * math.pi / 180.0, 0.3])
    # mu.append([-3.0 * math.pi / 180.0, 0.3])
    # mu.append([-2.0 * math.pi / 180.0, 0.3])
    # mu.append([-1.0 * math.pi / 180.0, 0.3])
    # mu.append([ 0.0 * math.pi / 180.0, 0.3])
    mu.append([ 1.0 * math.pi / 180.0, 0.3])
    for i in range(number_of_values):
         mu.append([values[i], 0.3])

    if plot_values:
        for i in range(len(values) + 2):
            plt.plot(np.round(mu[i][0] * 180 / math.pi,2) + 5, 0.3, 'bs')
        plt.xlabel('Alpha')
        plt.ylabel('Mach')
        plt.grid(True)
        plt.show()
    return mu

def get_multiple_params_mach():
    plot_values = True
    number_of_values = 4
    sampler = qmc.Halton(d=1)
    sample = sampler.random(number_of_values)
    #Mach infinit
    l_mach = [0.03]
    u_mach = [0.6]
    mu = []
    values = qmc.scale(sample, [l_mach[0]], [u_mach[0]])
    mu.append([0.0 * math.pi / 180.0, l_mach[0]])
    mu.append([0.0 * math.pi / 180.0, u_mach[0]])
    for i in range(number_of_values):
        mu.append([0.0 * math.pi / 180.0, values[i]])
    if plot_values:
        for i in range(len(values) + 2):
            plt.plot(mu[i][1], 0.0 * math.pi / 180.0 + 5.0, 'bs')
        plt.ylabel('Alpha')
        plt.xlabel('Mach')
        plt.grid(True)
        plt.show()
    return mu


def update_project_parameters(parameters, a_case,results_name):
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(a_case[0])
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(a_case[1])
    #storing results into a results folder
    #parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString(f"Results{results_name}[{str(np.round(a_case[1],2))}] - [{str(np.round(a_case[0]*180/math.pi,2)+5)}]")
    return parameters


def multiple_params_Train_Primal_ROM(mu):
    with open("ProjectParametersPrimalROM.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    parameters["output_processes"].RemoveValue("gid_output")
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
    for i in range(len(mu)):
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        print("::::::::::::Train HROM solution N° ",i,"::::::::::::")
        print("::::::::::::::::::::::::::::::::::::::::::::::::::")
        with open("ProjectParametersPrimalROM.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        if parameters["solver_settings"].Has("element_replace_settings"):
            parameters["solver_settings"].RemoveValue("element_replace_settings")
        parameters["output_processes"].RemoveValue("gid_output")
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
    parameters["output_processes"].RemoveValue("gid_output")
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


def primal_FOM(mach,angle):
    with open("ProjectParametersPrimalROM.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(angle)
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach)
    parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString(f"Results/FOM[{str(np.round(mach,2))}] - [{str(np.round(angle*180/math.pi,2)+5)}]")
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

def primal_ROM(mach,angle):
    with open("ProjectParametersPrimalROM.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(angle)
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach)
    parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString(f"Results/ROM[{str(np.round(mach,2))}] - [{str(np.round(angle*180/math.pi,2)+5)}]")
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

def primal_HROM(mach,angle):
    with open("ProjectParametersPrimalROM.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(angle)
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach)
    parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString(f"Results/HROM[{str(np.round(mach,2))}] - [{str(np.round(angle*180/math.pi,2)+5)}]")
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

def print_output_error(primal_model_part, output_parameters, primal_simulation, simulation):    
    gid_output = GiDOutputProcess(primal_model_part,'FluidModelPart', output_parameters)
    gid_output.ExecuteInitialize()
    gid_output.ExecuteBeforeSolutionLoop()
    for primal_node in primal_simulation._GetSolver().GetComputingModelPart().Nodes:
        primal_velocity = primal_node.GetSolutionStepValue(KratosMultiphysics.CompressiblePotentialFlowApplication.VELOCITY_POTENTIAL)
        primal_auxiliary_velocity = primal_node.GetSolutionStepValue(KratosMultiphysics.CompressiblePotentialFlowApplication.AUXILIARY_VELOCITY_POTENTIAL)
        rom_node = simulation._GetSolver().GetComputingModelPart().GetNode(primal_node.Id)
        rom_velocity = rom_node.GetSolutionStepValue(KratosMultiphysics.CompressiblePotentialFlowApplication.VELOCITY_POTENTIAL)
        rom_auxiliary_velocity = rom_node.GetSolutionStepValue(KratosMultiphysics.CompressiblePotentialFlowApplication.AUXILIARY_VELOCITY_POTENTIAL)
        if primal_velocity < 1e-10:
            primal_node.SetSolutionStepValue(KratosMultiphysics.CompressiblePotentialFlowApplication.VELOCITY_POTENTIAL, 0.0)
        else:
            primal_node.SetSolutionStepValue(KratosMultiphysics.CompressiblePotentialFlowApplication.VELOCITY_POTENTIAL, abs(primal_velocity-rom_velocity)/abs(primal_velocity)*100)
        if primal_auxiliary_velocity < 1e-10:
            primal_node.SetSolutionStepValue(KratosMultiphysics.CompressiblePotentialFlowApplication.AUXILIARY_VELOCITY_POTENTIAL, 0.0)
        else:
            primal_node.SetSolutionStepValue(KratosMultiphysics.CompressiblePotentialFlowApplication.AUXILIARY_VELOCITY_POTENTIAL, abs(primal_auxiliary_velocity-rom_auxiliary_velocity)/abs(primal_auxiliary_velocity)*100)
    gid_output.PrintOutput()

def run_test(angle, mach):
    fom_snapshots,tmfom = primal_FOM(mach,angle)
    setting_flags_rom_parameters(simulation_to_run = 'ROM', parameters_file_name = './PrimalRomParameters.json')
    rom_snapshots,tmrom = primal_ROM(mach,angle)
    setting_flags_rom_parameters(simulation_to_run = 'runHROM', parameters_file_name = './PrimalRomParameters.json')
    hrom_snapshots,tmhrom = primal_HROM(mach,angle)
    setting_flags_rom_parameters(simulation_to_run = 'ROM', parameters_file_name = './PrimalRomParameters.json')
    print("==========================> approximation error primal   FOM vs ROM: ",np.linalg.norm(fom_snapshots - rom_snapshots)/np.linalg.norm(fom_snapshots)*100,"%")
    print("==========================> approximation error primal  ROM vs HROM: ",np.linalg.norm(rom_snapshots - hrom_snapshots)/np.linalg.norm(rom_snapshots)*100,"%")
    print("time  FOM:", tmfom)
    print("time  ROM:", tmrom)
    print("time HROM:", tmhrom)

def error_calculate(angle, mach):
    with open("ProjectParametersPrimalROM.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(angle)
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach)
    
    parameters["output_processes"].RemoveValue("gid_output")
    #primal_parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString(f"Results/FOM[{str(np.round(mach,2))}] - [{str(np.round(angle*180/math.pi,2)+5)}]")
    
    primal_model = KratosMultiphysics.Model()
    primal_simulation = PotentialFlowAnalysis(primal_model, parameters)
    primal_simulation.Run()

    for process in primal_simulation._GetListOfOutputProcesses():
        if type(process) == KratosMultiphysics.RomApplication.calculate_rom_basis_output_process.CalculateRomBasisOutputProcess:
            primal_BasisOutputProcess = process
    primal_SnapshotsMatrix = primal_BasisOutputProcess._GetSnapshotsMatrix()


    with open("ProjectParametersPrimalROM.json",'r') as rom_parameter_file:
        rom_parameters = KratosMultiphysics.Parameters(rom_parameter_file.read())
    rom_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(angle)
    rom_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach)
    parameters["output_processes"].RemoveValue("gid_output")
    #primal_parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString(f"Results/ROM[{str(np.round(mach,2))}] - [{str(np.round(angle*180/math.pi,2)+5)}]")

    rom_parameters["output_processes"].RemoveValue("gid_output")
    if rom_parameters["solver_settings"].Has("element_replace_settings"):
        rom_parameters["solver_settings"].RemoveValue("element_replace_settings")
    model = KratosMultiphysics.Model()
    simulation = SetUpSimulationInstance(model,rom_parameters) 
    simulation.Run()

    for process in simulation._GetListOfOutputProcesses():
        if type(process) == KratosMultiphysics.RomApplication.calculate_rom_basis_output_process.CalculateRomBasisOutputProcess:
            BasisOutputProcess = process
    SnapshotsMatrix = BasisOutputProcess._GetSnapshotsMatrix()

    print("==========================> approximation error primal   FOM vs ROM: ",np.linalg.norm(primal_SnapshotsMatrix - SnapshotsMatrix)/np.linalg.norm(primal_SnapshotsMatrix)*100,"%")

    with open("ProjectParametersPrimal.json",'r') as output_parameter_file:
        output_parameters = KratosMultiphysics.Parameters(output_parameter_file.read())
    primal_model_part = primal_model.GetModelPart("MainModelPart")
    output_parameters["output_processes"].RemoveValue("rom_output")

    print_output_error(primal_model_part, output_parameters["output_processes"]["gid_output"][0]["Parameters"]["postprocess_parameters"], primal_simulation, simulation)



if __name__ == "__main__":

    create_basis = True

    if create_basis:
        KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Results')
        KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('ResultsFOM')
        KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('ResultsROM')
        KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('ResultsHROM')
        KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('PrimalRomParameters.json')
        KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('ROM test.post.lst')

        # mu = get_multiple_params()
        mu = get_multiple_params_angle()
        # mu = get_multiple_params_mach()
        
        
        primal_fom_snapshots = multiple_params_Train_Primal_ROM(mu)


        setting_flags_rom_parameters(simulation_to_run = 'trainHROM', parameters_file_name = './PrimalRomParameters.json')
        primal_rom_snapshots = multiple_params_Train_Primal_HROM(mu)
        setting_flags_rom_parameters(simulation_to_run = 'runHROM', parameters_file_name = './PrimalRomParameters.json')
        primal_hrom_snapshots = multiple_params_Primal_HROM(mu)
        setting_flags_rom_parameters(simulation_to_run = 'ROM', parameters_file_name = './PrimalRomParameters.json')
        print("==========================> approximation error primal   FOM vs ROM: ",np.linalg.norm(primal_fom_snapshots - primal_rom_snapshots)/np.linalg.norm(primal_fom_snapshots)*100,"%")
        print("==========================> approximation error primal  ROM vs HROM: ",np.linalg.norm(primal_rom_snapshots - primal_hrom_snapshots)/np.linalg.norm(primal_rom_snapshots)*100,"%")

    if not create_basis:
        # angle = -5.0 * math.pi / 180
        # mach  = 0.3
        # error_calculate(angle, mach)

        angle = -2.0 * math.pi / 180
        mach  = 0.3
        run_test(angle, mach)

