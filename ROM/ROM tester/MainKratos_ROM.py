import KratosMultiphysics
from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis import PotentialFlowAnalysis
from KratosMultiphysics.RomApplication.rom_testing_utilities import SetUpSimulationInstance
import KratosMultiphysics.kratos_utilities
import numpy as np
import math
import time


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



if __name__ == "__main__":
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Results')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('ResultsFOM')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('ResultsROM')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('ResultsHROM')

    angle = -0.5 * math.pi / 180
    mach  = 0.3

    fom_snapshots,tmfom = primal_FOM(mach,angle)

    rom_snapshots,tmrom = primal_ROM(mach,angle)

    hrom_snapshots,tmhrom = primal_HROM(mach,angle)

    print("==========================> approximation error primal   FOM vs ROM: ",np.linalg.norm(fom_snapshots - rom_snapshots)/np.linalg.norm(fom_snapshots)*100,"%")
    print("==========================> approximation error primal  ROM vs HROM: ",np.linalg.norm(rom_snapshots - hrom_snapshots)/np.linalg.norm(rom_snapshots)*100,"%")
    print("time  FOM:", tmfom)
    print("time  ROM:", tmrom)
    print("time HROM:", tmhrom)