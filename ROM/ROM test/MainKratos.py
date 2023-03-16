import KratosMultiphysics
from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis import PotentialFlowAnalysis
from KratosMultiphysics.RomApplication.rom_testing_utilities import SetUpSimulationInstance
import KratosMultiphysics.kratos_utilities
import numpy as np
from scipy.stats import qmc
import matplotlib.pyplot as plt
import compute_Cp_and_Cl
import math
import time
import json

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
    #parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString('Results/FOM'+str(mach))
    model = KratosMultiphysics.Model()
    simulation = PotentialFlowAnalysis(model, parameters)
    start=time.time()
    simulation.Run()

    plt.figure()
    plt.subplot(211)
    plt.title('FOM Naca 0012 Cp Distribution')
    nodal_Cp = np.zeros((len(model.GetModelPart("MainModelPart").Nodes)))
    for cond in model.GetModelPart("MainModelPart.Body2D_Body").Conditions:
        Cp = cond.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
        for node in cond.GetNodes():
            Cp_contribution = Cp/2.0
            nodal_Cp[node.Id] += Cp_contribution        
    for node in model.GetModelPart("MainModelPart.Body2D_Body").Nodes:
        plt.plot(node.X, nodal_Cp[node.Id], 'bs')
    plt.axis([-0.75, 0.75, -3.0, 1.5])
    plt.ylabel('Cp')
    plt.grid(True)
    plt.subplot(212)
    for node in model.GetModelPart("MainModelPart.Body2D_Body").Nodes:
        plt.plot(node.X, node.Y, 'ro')
    plt.axis([-0.75, 0.75, -0.75, 0.75])
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)
    plt.show(block=False)
    plt.savefig("FOM.png")

    # plt.title('Naca 0012 Cp Distribution')
    # nodal_Cp = np.zeros((len(model.GetModelPart("MainModelPart").Nodes)))
    # for cond in model.GetModelPart("MainModelPart.Body2D_Body").Conditions:
    #     Cp = cond.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
    #     for node in cond.GetNodes():
    #         Cp_contribution = Cp/2.0
    #         nodal_Cp[node.Id] += Cp_contribution  
    # for node in model.GetModelPart("MainModelPart.Body2D_Body").Nodes:
    #     plt.plot(node.X, nodal_Cp[node.Id], 'bs',node.X, node.Y, 'ro')
    # plt.axis([-0.75, 0.75, -2.5, 1.5])
    # plt.xlabel('x')
    # plt.ylabel('y  - Cp')
    # plt.grid(True)
    # plt.show()

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
    #parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString('Results/ROM'+str(mach))
    model = KratosMultiphysics.Model()
    simulation = SetUpSimulationInstance(model,parameters) 
    start=time.time()
    simulation.Run()

    plt.figure()
    plt.subplot(211)
    plt.title('ROM Naca 0012 Cp Distribution')
    nodal_Cp = np.zeros((len(model.GetModelPart("ModelPart").Nodes)))
    for cond in model.GetModelPart("MainModelPart.Body2D_Body").Conditions:
        Cp = cond.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
        for node in cond.GetNodes():
            Cp_contribution = Cp/2.0
            nodal_Cp[node.Id] += Cp_contribution        
    for node in model.GetModelPart("MainModelPart.Body2D_Body").Nodes:
        plt.plot(node.X, nodal_Cp[node.Id], 'bs')
    plt.axis([-0.75, 0.75, -3.0, 1.5])
    plt.ylabel('Cp')
    plt.grid(True)
    plt.subplot(212)
    for node in model.GetModelPart("MainModelPart.Body2D_Body").Nodes:
        plt.plot(node.X, node.Y, 'ro')
    plt.axis([-0.75, 0.75, -0.75, 0.75])
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)
    plt.show(block=False)
    plt.savefig("ROM.png")

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
    #parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString('Results/HROM'+str(mach))
    model = KratosMultiphysics.Model()
    simulation = SetUpSimulationInstance(model,parameters) 
    start=time.time()
    simulation.Run()

    plt.figure()
    plt.subplot(211)
    plt.title('HROM Naca 0012 Cp Distribution')
    nodal_Cp = np.zeros((len(model.GetModelPart("MainModelPart").Nodes)))
    for cond in model.GetModelPart("MainModelPart.Body2D_Body").Conditions:
        Cp = cond.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
        for node in cond.GetNodes():
            Cp_contribution = Cp/2.0
            nodal_Cp[node.Id] += Cp_contribution        
    for node in model.GetModelPart("MainModelPart.Body2D_Body").Nodes:
        plt.plot(node.X, nodal_Cp[node.Id], 'bs')
    plt.axis([-0.75, 0.75, -3.0, 1.5])
    plt.ylabel('Cp')
    plt.grid(True)
    plt.subplot(212)
    for node in model.GetModelPart("MainModelPart.Body2D_Body").Nodes:
        plt.plot(node.X, node.Y, 'ro')
    plt.axis([-0.75, 0.75, -0.75, 0.75])
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)
    plt.show(block=False)
    plt.savefig("HROM.png")

    end=time.time()  
    for process in simulation._GetListOfOutputProcesses():
        if type(process) == KratosMultiphysics.RomApplication.calculate_rom_basis_output_process.CalculateRomBasisOutputProcess:
            BasisOutputProcess = process
    SnapshotsMatrix = BasisOutputProcess._GetSnapshotsMatrix()
    tm = end - start
    return SnapshotsMatrix,tm




if __name__ == "__main__":
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('ROM test.post.lst')

    mach = 0.3

    fom_snapshots,tmfom = primal_FOM(mach)
    setting_flags_rom_parameters(simulation_to_run = 'ROM', parameters_file_name = './PrimalRomParameters.json')
    rom_snapshots,tmrom = primal_ROM(mach)
    setting_flags_rom_parameters(simulation_to_run = 'runHROM', parameters_file_name = './PrimalRomParameters.json')
    hrom_snapshots,tmhrom = primal_HROM(mach)
    setting_flags_rom_parameters(simulation_to_run = 'ROM', parameters_file_name = './PrimalRomParameters.json')
    print("==========================> approximation error primal   FOM vs ROM: ",np.linalg.norm(fom_snapshots - rom_snapshots)/np.linalg.norm(fom_snapshots)*100,"%")
    print("==========================> approximation error primal  ROM vs HROM: ",np.linalg.norm(rom_snapshots - hrom_snapshots)/np.linalg.norm(rom_snapshots)*100,"%")
    print("time  FOM:", tmfom)
    print("time  ROM:", tmrom)
    print("time HROM:", tmhrom)