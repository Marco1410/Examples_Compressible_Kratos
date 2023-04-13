import KratosMultiphysics
from KratosMultiphysics.RomApplication.rom_testing_utilities import SetUpSimulationInstance
from KratosMultiphysics.MeshMovingApplication.mesh_moving_analysis import MeshMovingAnalysis
from KratosMultiphysics.RomApplication.calculate_rom_basis_output_process import CalculateRomBasisOutputProcess
import KratosMultiphysics.kratos_utilities
import numpy as np
import importlib
import math
import time
import json

def CustomizeSimulation(cls, global_model, parameters):

    class CustomSimulation(cls):

        def __init__(self, model,project_parameters, custom_param = None):
            super().__init__(model,project_parameters)
            self.custom_param  = custom_param

        def Initialize(self):

            angle_of_attack = parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].GetDouble()
            parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(0.0) 
            with open("ProjectParametersMeshMoving.json",'r') as parameter_file:
                mesh_parameters = KratosMultiphysics.Parameters(parameter_file.read()) 
            mesh_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["rotation_angle"].SetDouble(angle_of_attack)
            mesh_simulation = MeshMovingAnalysis(self.model,mesh_parameters)

            super().Initialize()
            
            mesh_simulation.Run()

        def FinalizeSolutionStep(self):
            super().FinalizeSolutionStep()

        def CustomMethod(self):
            return self.custom_param

    return CustomSimulation(global_model, parameters)

def _GetAnalysisStageClass(parameters):
    analysis_stage_module_name = parameters["analysis_stage"].GetString()
    analysis_stage_class_name = analysis_stage_module_name.split('.')[-1]
    analysis_stage_class_name = ''.join(x.title() for x in analysis_stage_class_name.split('_'))
    analysis_stage_module = importlib.import_module(analysis_stage_module_name)
    analysis_stage_class = getattr(analysis_stage_module, analysis_stage_class_name)
    return analysis_stage_class

def _ChangeRomFlags(simulation_to_run = 'ROM'):
        parameters_file_name = './RomParameters.json'
        with open(parameters_file_name, 'r+') as parameter_file:
            f=json.load(parameter_file)
            if simulation_to_run=='ROM':
                f['train_hrom']=False
                f['run_hrom']=False
            elif simulation_to_run=='HROM':
                f['train_hrom']=False
                f['run_hrom']=True
            else:
                raise Exception(f'Unknown flag "{simulation_to_run}" change for RomParameters.json')
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
       
    start=time.time()
    analysis_stage_class = _GetAnalysisStageClass(parameters)
    simulation = CustomizeSimulation(analysis_stage_class,model,parameters)
    simulation.Run()
    end=time.time()  

    SnapshotsMatrix = []
    for process in simulation._GetListOfOutputProcesses():
        if isinstance(process, CalculateRomBasisOutputProcess):
            BasisOutputProcess = process
    SnapshotsMatrix.append(BasisOutputProcess._GetSnapshotsMatrix()) 
    SnapshotsMatrix = np.block(SnapshotsMatrix)

    tm = end - start

    return SnapshotsMatrix,tm



def primal_ROM(mach,angle):
    with open("ProjectParametersPrimalROM.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(angle)
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach)
    parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString(f"Results/ROM[{str(np.round(mach,2))}] - [{str(np.round(angle*180/math.pi,2)+5)}]")
    _ChangeRomFlags(simulation_to_run = "ROM")
    model = KratosMultiphysics.Model()
        
    start=time.time()
    analysis_stage_class = type(SetUpSimulationInstance(model, parameters))
    simulation = CustomizeSimulation(analysis_stage_class,model,parameters)
    simulation.Run()
    end=time.time()  
    
    # SnapshotsMatrix = []
    # for process in simulation._GetListOfOutputProcesses():
    #     if isinstance(process, CalculateRomBasisOutputProcess):
    #         BasisOutputProcess = process
    # SnapshotsMatrix.append(BasisOutputProcess._GetSnapshotsMatrix()) 
    # SnapshotsMatrix = np.block(SnapshotsMatrix)

    tm = end - start

    return tm



def primal_HROM(mach,angle):
    with open("ProjectParametersPrimalROM.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(angle)
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach)
    parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString(f"Results/HROM[{str(np.round(mach,2))}] - [{str(np.round(angle*180/math.pi,2)+5)}]")
    _ChangeRomFlags(simulation_to_run = "HROM")
    model = KratosMultiphysics.Model()
        
    start=time.time()
    analysis_stage_class = type(SetUpSimulationInstance(model, parameters))
    simulation = CustomizeSimulation(analysis_stage_class,model,parameters)
    simulation.Run()
    end=time.time()  
    
    # SnapshotsMatrix = []
    # for process in simulation._GetListOfOutputProcesses():
    #     if isinstance(process, CalculateRomBasisOutputProcess):
    #         BasisOutputProcess = process
    # SnapshotsMatrix.append(BasisOutputProcess._GetSnapshotsMatrix()) 
    # SnapshotsMatrix = np.block(SnapshotsMatrix)

    tm = end - start

    return tm



if __name__ == "__main__":
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Results')

    angle = 0.5 * math.pi / 180
    mach  = 0.3

    tmfom  = primal_FOM(mach,angle)
    tmrom  = primal_ROM(mach,angle)
    tmhrom = primal_HROM(mach,angle)

    #print("==========================> approximation error primal   FOM vs ROM: ",np.linalg.norm(fom_snapshots - rom_snapshots)/np.linalg.norm(fom_snapshots)*100,"%")
    # print("==========================> approximation error primal  ROM vs HROM: ",np.linalg.norm(rom_snapshots - hrom_snapshots)/np.linalg.norm(rom_snapshots)*100,"%")
    print("time  FOM:", tmfom)
    print("time  ROM:", tmrom)
    print("time HROM:", tmhrom)