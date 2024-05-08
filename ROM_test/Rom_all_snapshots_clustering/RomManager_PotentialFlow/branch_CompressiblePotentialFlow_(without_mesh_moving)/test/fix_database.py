import os
import sys
import time
import importlib
import numpy as np
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

import KratosMultiphysics

def CreateAnalysisStageWithFlushInstance(cls, global_model, parameters):
    class AnalysisStageWithFlush(cls):

        def __init__(self, model,project_parameters, flush_frequency=10.0):
            super().__init__(model,project_parameters)
            self.flush_frequency = flush_frequency
            self.last_flush = time.time()
            sys.stdout.flush()

        def InitializeSolutionStep(self):
            super().InitializeSolutionStep()
            sys.stdout.flush()

        def FinalizeSolutionStep(self):
            super().FinalizeSolutionStep()
    
            if self.parallel_type == "OpenMP":
                now = time.time()
                if now - self.last_flush > self.flush_frequency:
                    sys.stdout.flush()
                    self.last_flush = now

    return AnalysisStageWithFlush(global_model, parameters)

if __name__ == "__main__":

    dir = "DataBase/"
    fix_database = False
    plot_values  = False

    with os.scandir(dir) as files:
        files = [file.name for file in files if file.is_file() and file.name.endswith('.npy')]

    print(" ")
    print(len(files)," files found")

    if plot_values:
        data_set_m = np.zeros(len(files))
        data_set_a = np.zeros(len(files))
        fig = plt.figure()
        fig.set_figwidth(12.0)
        fig.set_figheight(8.0)

    for n,file in enumerate(files):
        name = file.removesuffix(".npy")
        angle_of_attack = np.double(name.split(",")[0])
        mach_infinity   = np.double(name.split(",")[1].removeprefix(" "))

        if plot_values:
            data_set_a[n] = angle_of_attack
            data_set_m[n] = mach_infinity

        if not os.path.exists(dir+"/Data/"+name+".dat") or not os.path.exists(dir+"/Data/full_"+name+".dat") or not os.path.exists(dir+"/Gid_Results/"+name+".post.bin"):

            print(" ")
            print("case: ", name)

            if fix_database:
                with open("ProjectParameters.json", 'r') as parameter_file:
                    parameters = KratosMultiphysics.Parameters(parameter_file.read())

                parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(angle_of_attack)
                parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach_infinity)
                parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString("DataBase/Gid_Results/"+ name)

                analysis_stage_module_name = parameters["analysis_stage"].GetString()
                analysis_stage_class_name = analysis_stage_module_name.split('.')[-1]
                analysis_stage_class_name = ''.join(x.title() for x in analysis_stage_class_name.split('_'))

                analysis_stage_module = importlib.import_module(analysis_stage_module_name)
                analysis_stage_class = getattr(analysis_stage_module, analysis_stage_class_name)
                
                global_model = KratosMultiphysics.Model()
                simulation = CreateAnalysisStageWithFlushInstance(analysis_stage_class, global_model, parameters)
                simulation.Run()

                skin_data = "DataBase/Data/" + name + ".dat"
                fout = open(skin_data,'w')

                modelpart = global_model["MainModelPart.Body2D_Body"]
                for node in modelpart.Nodes:
                    x = node.X ; y = node.Y ; z = node.Z
                    cp = node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
                    fout.write("%s %s %s %s\n" %(x,y,z,cp))
                fout.close()

                fout=open("DataBase/Data/full_" + name + ".dat",'w')
                modelpart = global_model["MainModelPart"]
                for node in modelpart.Nodes:
                    id_node = node.Id
                    velocity_potential = node.GetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL)
                    auxiliary_velocity_potential = node.GetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL)
                    fout.write("%s %s %s\n" %(id_node, velocity_potential, auxiliary_velocity_potential))
                fout.close()

    if plot_values:
        fig = plt.plot(data_set_m, data_set_a, "ob", label="Data Base")
        fig = plt.title('Mach vs Angle of attack')
        fig = plt.ylabel('Alpha')
        fig = plt.xlabel('Mach')
        fig = plt.grid(True)
        fig = plt.legend(bbox_to_anchor=(.85, 1.03, 1., .102), loc='upper left', borderaxespad=0.)
        fig = plt.savefig(dir+"0.0, 0.0-DataBase.png")
        fig = plt.close('all')

        print(" ")
        print("DB plot ready")
        print(" ")