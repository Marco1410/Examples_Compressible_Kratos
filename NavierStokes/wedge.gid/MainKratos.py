import sys
import time

import KratosMultiphysics
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

class FluidDynamicsAnalysisWithFlush(FluidDynamicsAnalysis):

    def __init__(self, model, project_parameters, flush_frequency=10.0):
        super().__init__(model, project_parameters)
        self.flush_frequency = flush_frequency
        self.last_flush = time.time()
        sys.stdout.flush()

    def Initialize(self):
        super().Initialize()
        sys.stdout.flush()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        if self.parallel_type == "OpenMP":
            now = time.time()
            if now - self.last_flush > self.flush_frequency:
                sys.stdout.flush()
                self.last_flush = now
    
    # def _GetOrderOfProcessesInitialization(self):
    #     """This function is overridden in order to set 
    #     the initialization order of the processes.
    #     """
    #     return ["initial_conditions_process_list", "boundary_conditions_process_list", "mesh_adaptivity_processes"]

if __name__ == "__main__":

    with open("ProjectParameters.json", 'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    global_model = KratosMultiphysics.Model()
    simulation = FluidDynamicsAnalysisWithFlush(global_model, parameters)
    simulation.Run()
