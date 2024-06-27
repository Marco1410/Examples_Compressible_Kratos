import sys
import time
import importlib

import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp 

def CreateAnalysisStageWithFlushInstance(cls, global_model, parameters):
    class AnalysisStageWithFlush(cls):

        def __init__(self, model,project_parameters, flush_frequency=10.0):
            super().__init__(model,project_parameters)
            self.flush_frequency = flush_frequency
            self.last_flush = time.time()
            sys.stdout.flush()
        
        @classmethod
        def distance_to_corner_top(cls, u: KratosMultiphysics.Node):
            return (u.X - 36.467168)**2 + (u.Y - 3.006712)**2 + (u.Z - 3.513132)**2

        def _find_root_node_top(self):
            eps_2 = 1e-6
            for node in self._GetSolver().GetComputingModelPart().Nodes:
                if self.distance_to_corner_top(node) < eps_2:
                    self._root_node_top = node.Id
                    return
            raise RuntimeError("Root node top not found")
        
        @classmethod
        def distance_to_corner_bottom(cls, u: KratosMultiphysics.Node):
            return (u.X - 36.466426)**2 + (u.Y - 3.006438)**2 + (u.Z - 3.500111)**2

        def _find_root_node_bottom(self):
            eps_2 = 1e-6
            for node in self._GetSolver().GetComputingModelPart().Nodes:
                if self.distance_to_corner_bottom(node) < eps_2:
                    self._root_node_bottom = node.Id
                    return
            raise RuntimeError("Root node bottom not found")
            
        def Initialize(self):
            super().Initialize()
            sys.stdout.flush()
            # self._find_root_node_top()
            # self._find_root_node_bottom()
        
        def InitializeSolutionStep(self):
            # node = self._GetSolver().GetComputingModelPart().GetNode(self._root_node_top)
            # node.SetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL, 0.0)
            # node.Fix(CPFApp.VELOCITY_POTENTIAL)
            # node = self._GetSolver().GetComputingModelPart().GetNode(self._root_node_bottom)
            # node.SetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL, 0.0)
            # node.Fix(CPFApp.VELOCITY_POTENTIAL)
            # node.SetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL, 200.0)
            # node.Fix(CPFApp.AUXILIARY_VELOCITY_POTENTIAL)
            super().InitializeSolutionStep()

        def FinalizeSolutionStep(self):
            super().FinalizeSolutionStep()

            if self.parallel_type == "OpenMP":
                now = time.time()
                if now - self.last_flush > self.flush_frequency:
                    sys.stdout.flush()
                    self.last_flush = now

    return AnalysisStageWithFlush(global_model, parameters)

if __name__ == "__main__":

    with open("ProjectParameters.json", 'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    analysis_stage_module_name = parameters["analysis_stage"].GetString()
    analysis_stage_class_name = analysis_stage_module_name.split('.')[-1]
    analysis_stage_class_name = ''.join(x.title() for x in analysis_stage_class_name.split('_'))

    analysis_stage_module = importlib.import_module(analysis_stage_module_name)
    analysis_stage_class = getattr(analysis_stage_module, analysis_stage_class_name)

    global_model = KratosMultiphysics.Model()
    simulation = CreateAnalysisStageWithFlushInstance(analysis_stage_class, global_model, parameters)
    simulation.Run()
