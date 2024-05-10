import os
import sys
import time
import importlib
import numpy as np
from scipy.stats import qmc
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis import PotentialFlowAnalysis
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

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# get multiple parameters
#

def get_multiple_params_by_Halton_sequence(number_of_values,angle,mach,fix_corners_of_parametric_space):
    if np.abs(angle[1]-angle[0])< 1e-2:
        if fix_corners_of_parametric_space and number_of_values < 2:
            print("Setting number of values to 2.")
            number_of_values = 2
        sampler = qmc.Halton(d=1)
        # sampler = qmc.LatinHypercube(d=1)
        mu = []
        if number_of_values > 0:
            sample = sampler.random(number_of_values)
            values = qmc.scale(sample, [mach[0]],[mach[1]])
            if fix_corners_of_parametric_space and number_of_values >= 2:
                values[0] = mach[0]
                values[number_of_values-1] = mach[1]
            for i in range(number_of_values):
                #Angle of attack , Mach infinit
                mu.append([angle[0], values[i]])
    elif np.abs(mach[1]-mach[0])< 1e-3:
        if fix_corners_of_parametric_space and number_of_values < 2:
            print("Setting number of values to 2.")
            number_of_values = 2
        sampler = qmc.Halton(d=1)
        # sampler = qmc.LatinHypercube(d=1)
        mu = []
        if number_of_values > 0:
            sample = sampler.random(number_of_values)
            values = qmc.scale(sample, [angle[0]],[angle[1]])
            if fix_corners_of_parametric_space and number_of_values >= 2:
                values[0] = angle[0]
                values[number_of_values-1] = angle[1]
            for i in range(number_of_values):
                #Angle of attack , Mach infinit
                mu.append([values[i], mach[0]])
    else:
        if fix_corners_of_parametric_space and number_of_values < 4:
            print("Setting number of values to 4.")
            number_of_values = 4
        sampler = qmc.Halton(d=2)
        # sampler = qmc.LatinHypercube(d=2)
        mu = []
        if number_of_values > 0:
            sample = sampler.random(number_of_values)
            values = qmc.scale(sample, [angle[0],mach[0]], [angle[1],mach[1]])
            if fix_corners_of_parametric_space and number_of_values >= 4:
                values[0,0] = angle[0]
                values[0,1] = mach[0]
                values[1,0] = angle[0]
                values[1,1] = mach[1]
                values[number_of_values-1,0] = angle[1]
                values[number_of_values-1,1] = mach[1]
                values[number_of_values-2,0] = angle[1]
                values[number_of_values-2,1] = mach[0]
            for i in range(number_of_values):
                #Angle of attack , Mach infinit
                mu.append([values[i,0], values[i,1]])
    return mu

if __name__ == "__main__":

    dir = "DataBase/"
    create       = False
    fix          = False
    plot_values  = False
    
    NumberofPoints = 50
    mach_range     = [ 0.72, 0.73]
    angle_range    = [ 1.00, 2.00]

    if create:
        fix_corners_of_parametric_space = False
        mu_values = get_multiple_params_by_Halton_sequence(NumberofPoints, angle_range, mach_range,
                                                           fix_corners_of_parametric_space)
        for n,mu in enumerate(mu_values):
            angle_of_attack = mu[0]
            mach_infinity   = mu[1]
            name = str(angle_of_attack) + ", " + str(mach_infinity)

            if not os.path.exists(dir+"/"+name+".npy"):
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

                skin_data = dir+"/Data/" + name + ".dat"
                fout = open(skin_data,'w')
                modelpart = global_model["MainModelPart.Body2D_Body"]
                for node in modelpart.Nodes:
                    x = node.X ; y = node.Y ; z = node.Z
                    cp = node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
                    fout.write("%s %s %s %s\n" %(x,y,z,cp))
                fout.close()
                full_data_name = dir+"/Data/full_" + name + ".dat"
                fout=open(full_data_name,'w')
                modelpart = global_model["MainModelPart"]
                for node in modelpart.Nodes:
                    id_node = node.Id
                    velocity_potential = node.GetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL)
                    auxiliary_velocity_potential = node.GetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL)
                    fout.write("%s %s %s\n" %(id_node, velocity_potential, auxiliary_velocity_potential))
                fout.close()

                skin_data_uncorrected_case_name = dir+"/NotSharpedSolution/SkinData/"+name+".dat"
                parameters["solver_settings"]["scheme_settings"]["update_critical_mach"].SetDouble(0.85)
                parameters["solver_settings"]["scheme_settings"]["update_upwind_factor_constant"].SetDouble(2.0)
                parameters["solver_settings"]["scheme_settings"]["update_transonic_tolerance"].SetDouble(1e-30)
                if parameters["output_processes"].Has("gid_output"):
                    parameters["output_processes"].RemoveValue("gid_output")
                aux_model = KratosMultiphysics.Model()
                aux_simulation = PotentialFlowAnalysis(aux_model,parameters)
                aux_simulation.Run()
                fout = open(skin_data_uncorrected_case_name,'w')
                modelpart = aux_model["MainModelPart.Body2D_Body"]
                for node in modelpart.Nodes:
                    x = node.X ; y = node.Y ; z = node.Z
                    cp = node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
                    fout.write("%s %s %s %s\n" %(x,y,z,cp))
                fout.close()

    with os.scandir(dir) as files:
        files = [file.name for file in files if file.is_file() and file.name.endswith('.npy')]

    print(" ")
    print(len(files)," files found")
    count = 0

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

        if (   not os.path.exists(dir+"/Data/"+name+".dat") 
            or not os.path.exists(dir+"/Data/full_"+name+".dat") 
            or not os.path.exists(dir+"/Gid_Results/"+name+".post.bin")
            or not os.path.exists(dir+"/NotSharpedSolution/SkinData/"+name+".dat")):
            count += 1
            if fix:
                print(" ")
                print("case: ", name)

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

                skin_data = dir+"/Data/" + name + ".dat"
                fout = open(skin_data,'w')

                modelpart = global_model["MainModelPart.Body2D_Body"]
                for node in modelpart.Nodes:
                    x = node.X ; y = node.Y ; z = node.Z
                    cp = node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
                    fout.write("%s %s %s %s\n" %(x,y,z,cp))
                fout.close()

                fout=open(dir+"/Data/full_" + name + ".dat",'w')
                modelpart = global_model["MainModelPart"]
                for node in modelpart.Nodes:
                    id_node = node.Id
                    velocity_potential = node.GetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL)
                    auxiliary_velocity_potential = node.GetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL)
                    fout.write("%s %s %s\n" %(id_node, velocity_potential, auxiliary_velocity_potential))
                fout.close()

                skin_data_name = dir+"/NotSharpedSolution/SkinData/"+name+".dat"
                if not os.path.exists(skin_data_name):
                    parameters["solver_settings"]["scheme_settings"]["update_critical_mach"].SetDouble(0.85)
                    parameters["solver_settings"]["scheme_settings"]["update_upwind_factor_constant"].SetDouble(2.0)
                    parameters["solver_settings"]["scheme_settings"]["update_transonic_tolerance"].SetDouble(1e-30)
                    if parameters["output_processes"].Has("gid_output"):
                        parameters["output_processes"].RemoveValue("gid_output")
                    aux_model = KratosMultiphysics.Model()
                    aux_simulation = PotentialFlowAnalysis(aux_model,parameters)
                    aux_simulation.Run()
                    fout = open(skin_data_name,'w')
                    modelpart = aux_model["MainModelPart.Body2D_Body"]
                    for node in modelpart.Nodes:
                        x = node.X ; y = node.Y ; z = node.Z
                        cp = node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
                        fout.write("%s %s %s %s\n" %(x,y,z,cp))
                    fout.close()

    print(count, " file/s to fix.")

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