import sys
import time
import math

import KratosMultiphysics
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

class FluidDynamicsAnalysisWithFlush(FluidDynamicsAnalysis):

    def __init__(self, model, project_parameters, flush_frequency=10.0):
        super().__init__(model, project_parameters)
        self.flush_frequency = flush_frequency
        self.last_flush = time.time()
        sys.stdout.flush()
        
    @classmethod
    def distance_to_corner(cls, u: KratosMultiphysics.Node):
        return (u.X - 0.6)**2 + (u.Y - 0.0)**2 + (u.Y - 0.0)**2

    def _find_lower_corner_node(self):
        eps_2 = 1e-16
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            if self.distance_to_corner(node) < eps_2:
                self._lower_corner_node = node.Id
                return
        raise RuntimeError("Lower corner node not found")

    def Initialize(self):
        super().Initialize()
        sys.stdout.flush()
        
    def InitializeSolutionStep(self):
        self._find_lower_corner_node()
        node = self._GetSolver().GetComputingModelPart().GetNode(self._lower_corner_node)
        node.SetSolutionStepValue(KratosMultiphysics.MOMENTUM_X, 0.0)
        node.SetSolutionStepValue(KratosMultiphysics.MOMENTUM_Y, 0.0)
        node.Fix(KratosMultiphysics.MOMENTUM_X)
        node.Fix(KratosMultiphysics.MOMENTUM_Y)
        super().InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        self.Flagg()
    
    def Flagg(self):
        step = self._GetSolver().GetComputingModelPart().ProcessInfo.GetValue(KratosMultiphysics.STEP)
        if  step > 1 :
            e_totalEnergy = 0
            n_totalEnergy = 0
            for node in self._GetSolver().GetComputingModelPart().Nodes:
                e_totalEnergy = e_totalEnergy + (node.GetSolutionStepValue(KratosMultiphysics.TOTAL_ENERGY,1)-node.GetSolutionStepValue(KratosMultiphysics.TOTAL_ENERGY))**2
                n_totalEnergy = n_totalEnergy + (node.GetSolutionStepValue(KratosMultiphysics.TOTAL_ENERGY,1))**2 

            e_totalEnergy = (e_totalEnergy/n_totalEnergy)**0.5
                
            if math.isnan(e_totalEnergy):
                KratosMultiphysics.Logger.Print("               NaN value in TOTAL ENERGY EQUATION")
                sys.exit('               STOP')

        if self.parallel_type == "OpenMP":
            now = time.time()
            if now - self.last_flush > self.flush_frequency:
                sys.stdout.flush()
                self.last_flush = now
    
    def _GetOrderOfProcessesInitialization(self):
        """This function is overridden in order to set 
        the initialization order of the processes.
        """
        return ["mesh_adaptivity_processes", "initial_conditions_process_list", "boundary_conditions_process_list"]

if __name__ == "__main__":

    with open("ProjectParameters.json", 'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    global_model = KratosMultiphysics.Model()
    simulation = FluidDynamicsAnalysisWithFlush(global_model, parameters)
    simulation.Run()
