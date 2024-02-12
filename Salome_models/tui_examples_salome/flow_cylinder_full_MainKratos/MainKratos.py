import KratosMultiphysics
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

import os
import sys
import time

class FluidDynamicsAnalysisWithFlush(FluidDynamicsAnalysis):

    def __init__(self,model,project_parameters,flush_frequency=10.0):
        super(FluidDynamicsAnalysisWithFlush,self).__init__(model,project_parameters)
        self.flush_frequency = flush_frequency
        self.last_flush = time.time()

    def FinalizeSolutionStep(self):
        super(FluidDynamicsAnalysisWithFlush,self).FinalizeSolutionStep()

        if self.parallel_type == "OpenMP":
            now = time.time()
            if now - self.last_flush > self.flush_frequency:
                sys.stdout.flush()
                self.last_flush = now

def create_salome_mesh(salome_script_name):
    import time as time
    import subprocess

    start_time = time.time()

    salome_cmd = "salome -t python"
    salome_exe = " ".join([salome_cmd, salome_script_name])
    sp = subprocess.Popen(["/bin/bash", "-i", "-c", salome_exe])
    sp.communicate()

    exe_time = time.time() - start_time

    print('Executing SALOME with "' + salome_script_name +
        '" took ' + str(round(exe_time, 2)) + ' sec')

if __name__ == "__main__":

    if not os.path.exists("mesh.med"):
        create_salome_mesh("salome_model.py")

    with open("ProjectParameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = FluidDynamicsAnalysisWithFlush(model,parameters)
    simulation.Run()
