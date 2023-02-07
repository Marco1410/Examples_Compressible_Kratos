import sys
import time
import importlib

import KratosMultiphysics

def CreateAnalysisStageWithFlushInstance(cls, global_model, parameters):
    class AnalysisStageWithFlush(cls):

        def __init__(self, model,project_parameters, flush_frequency=10.0):
            super().__init__(model,project_parameters)
            self.flush_frequency = flush_frequency
            self.last_flush = time.time()
            sys.stdout.flush()
         
        @classmethod
        def distance_to_corner1(cls, u: KratosMultiphysics.Node):
            return (u.X - 0.0)**2 + (u.Y - 1.0)**2

        def _find_corner_node1(self):
            eps_2 = 1e-16
            for node in self._GetSolver().GetComputingModelPart().Nodes:
                if self.distance_to_corner1(node) < eps_2:
                    self._corner_node1 = node.Id
                    KratosMultiphysics.Logger.Print(f"First corner node has ID: {self._corner_node1}")
                    return
            raise RuntimeError("First corner node not found")
                    
        @classmethod
        def distance_to_corner2(cls, u: KratosMultiphysics.Node):
            return (u.X - 1.0)**2 + (u.Y - 1.0)**2 

        def _find_corner_node2(self):
            eps_2 = 1e-16
            for node in self._GetSolver().GetComputingModelPart().Nodes:
                if self.distance_to_corner2(node) < eps_2:
                    self._corner_node2 = node.Id
                    KratosMultiphysics.Logger.Print(f"Second corner node has ID: {self._corner_node2}")
                    return
            raise RuntimeError("Second corner node not found")
                  
        @classmethod
        def distance_to_corner3(cls, u: KratosMultiphysics.Node):
            return (u.X - 1.0)**2 + (u.Y - 0.0)**2 

        def _find_corner_node3(self):
            eps_2 = 1e-16
            for node in self._GetSolver().GetComputingModelPart().Nodes:
                if self.distance_to_corner3(node) < eps_2:
                    self._corner_node3 = node.Id
                    KratosMultiphysics.Logger.Print(f"Third corner node has ID: {self._corner_node3}")
                    return
            raise RuntimeError("Third corner node not found")
                
        @classmethod
        def distance_to_corner4(cls, u: KratosMultiphysics.Node):
            return (u.X - 0.0)**2 + (u.Y - 0.0)**2 

        def _find_corner_node4(self):
            eps_2 = 1e-16
            for node in self._GetSolver().GetComputingModelPart().Nodes:
                if self.distance_to_corner4(node) < eps_2:
                    self._corner_node4 = node.Id
                    KratosMultiphysics.Logger.Print(f"Fourth corner node has ID: {self._corner_node4}")
                    return
            raise RuntimeError("Fourth corner node not found")
                   
        def Initialize(self):
            super().Initialize()
            sys.stdout.flush()
            self._find_corner_node1()
            self._find_corner_node2()
            self._find_corner_node3()
            self._find_corner_node4()
        
        def InitializeSolutionStep(self):
            node = self._GetSolver().GetComputingModelPart().GetNode(self._corner_node1)
            node.SetSolutionStepValue(KratosMultiphysics.MOMENTUM_X, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.MOMENTUM_Y, 0.0)
            node.Fix(KratosMultiphysics.MOMENTUM_X)
            node.Fix(KratosMultiphysics.MOMENTUM_Y)
            node = self._GetSolver().GetComputingModelPart().GetNode(self._corner_node2)
            node.SetSolutionStepValue(KratosMultiphysics.MOMENTUM_X, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.MOMENTUM_Y, 0.0)
            node.Fix(KratosMultiphysics.MOMENTUM_X)
            node.Fix(KratosMultiphysics.MOMENTUM_Y)
            node = self._GetSolver().GetComputingModelPart().GetNode(self._corner_node3)
            node.SetSolutionStepValue(KratosMultiphysics.MOMENTUM_X, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.MOMENTUM_Y, 0.0)
            node.Fix(KratosMultiphysics.MOMENTUM_X)
            node.Fix(KratosMultiphysics.MOMENTUM_Y)
            node = self._GetSolver().GetComputingModelPart().GetNode(self._corner_node4)
            node.SetSolutionStepValue(KratosMultiphysics.MOMENTUM_X, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.MOMENTUM_Y, 0.0)
            node.Fix(KratosMultiphysics.MOMENTUM_X)
            node.Fix(KratosMultiphysics.MOMENTUM_Y)
            super().InitializeSolutionStep()
            
        def FinalizeSolutionStep(self):
            super().FinalizeSolutionStep()
 
        def Flush(self):
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
