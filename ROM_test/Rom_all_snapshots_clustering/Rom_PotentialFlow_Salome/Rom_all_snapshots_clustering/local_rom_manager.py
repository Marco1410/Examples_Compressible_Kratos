import KratosMultiphysics
from KratosMultiphysics.RomApplication.rom_testing_utilities import SetUpSimulationInstance
from KratosMultiphysics.RomApplication.calculate_rom_basis_output_process import CalculateRomBasisOutputProcess
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition
import numpy as np
import os
import importlib
import json
from matplotlib import pyplot as plt
from pathlib import Path

#import kmeans
from sklearn.cluster import KMeans
######



class LocalRomManager(object):

    def __init__(self,project_parameters_name, general_rom_manager_parameters, CustomizeSimulation, UpdateProjectParameters):
        self.project_parameters_name = project_parameters_name
        self.general_rom_manager_parameters = general_rom_manager_parameters
        self.CustomizeSimulation = CustomizeSimulation
        self.UpdateProjectParameters = UpdateProjectParameters
        self.ROMvsFOM_train=self.ROMvsHROM_train=self.ROMvsFOM_test=self.ROMvsHROM_test=0.0
        self.train_info_steps_list = []
        self.test_info_steps_list = []


    def Fit(self, mu_train=[None], use_non_linear_steps=False, clustering_method="k-means", n_clusters=5, bisection_tolerance=0.15,  POD_tolerance=1e-12, use_aux_step=False):
        mu_train_errors = []; cluster_list = []
        if (len(mu_train) > 0):
            chosen_projection_strategy = self.general_rom_manager_parameters["projection_strategy"].GetString()
            training_stages = self.general_rom_manager_parameters["rom_stages_to_train"].GetStringArray()
            #######################
            ######  Galerkin ######
            if chosen_projection_strategy == "galerkin":
                if any(item == "ROM" for item in training_stages):
                    fom_snapshots, n_clusters = self.__LaunchTrainROM(mu_train, use_non_linear_steps, cluster_list, clustering_method, n_clusters, bisection_tolerance,  POD_tolerance, use_aux_step)
                    # self._ChangeRomFlags(simulation_to_run = "GalerkinROM")
                    rom_snapshots = self.__LaunchROM(mu_train, cluster_list, "GalerkinROM")
                    self.ROMvsFOM_train = np.linalg.norm(fom_snapshots - rom_snapshots)*100/ np.linalg.norm(fom_snapshots)
                    for i in range(len(mu_train)):
                        error = np.linalg.norm(fom_snapshots[:,i] - rom_snapshots[:,i])*100/ np.linalg.norm(fom_snapshots[:,i])
                        mu_train_errors.append([mu_train[i][0],mu_train[i][1],error,-1])

                if any(item == "HROM" for item in training_stages):
                    #FIXME there will be an error if we only train HROM, but not ROM
                    # self._ChangeRomFlags(simulation_to_run = "trainHROMGalerkin")
                    self.__LaunchTrainHROM(mu_train, n_clusters, cluster_list, "trainHROMGalerkin")
                    # self._ChangeRomFlags(simulation_to_run = "runHROMGalerkin")
                    hrom_snapshots = self.__LaunchHROM(mu_train, cluster_list, "runHROMGalerkin")
                    self.ROMvsHROM_train = np.linalg.norm(rom_snapshots - hrom_snapshots)*100/ np.linalg.norm(rom_snapshots)
                    for i in range(len(mu_train)):
                        error = np.linalg.norm(rom_snapshots[:,i] - hrom_snapshots[:,i])*100/ np.linalg.norm(rom_snapshots[:,i])
                        mu_train_errors[i][3] = error
            #######################

            #######################################
            ##  Least-Squares Petrov Galerkin   ###
            elif chosen_projection_strategy == "lspg":
                if any(item == "ROM" for item in training_stages):
                    fom_snapshots, n_clusters = self.__LaunchTrainROM(mu_train, use_non_linear_steps, cluster_list, clustering_method, n_clusters, bisection_tolerance,  POD_tolerance, use_aux_step)
                    # self._ChangeRomFlags(simulation_to_run = "lspg")
                    rom_snapshots = self.__LaunchROM(mu_train, cluster_list, "lspg")
                    self.ROMvsFOM_train = np.linalg.norm(fom_snapshots - rom_snapshots)*100/ np.linalg.norm(fom_snapshots)
                    for i in range(len(mu_train)):
                        error = np.linalg.norm(fom_snapshots[:,i] - rom_snapshots[:,i])*100/ np.linalg.norm(fom_snapshots[:,i])
                        mu_train_errors.append([mu_train[i][0],mu_train[i][1],error,-1])
                if any(item == "HROM" for item in training_stages):
                    # Change the flags to train the HROM for LSPG
                    # self._ChangeRomFlags(simulation_to_run = "trainHROMLSPG")
                    self.__LaunchTrainHROM(mu_train, n_clusters, cluster_list, "trainHROMLSPG")

                    # Change the flags to run the HROM for LSPG
                    # self._ChangeRomFlags(simulation_to_run = "runHROMLSPG")
                    hrom_snapshots = self.__LaunchHROM(mu_train, cluster_list, "runHROMLSPG")
                    self.ROMvsHROM_train = np.linalg.norm(rom_snapshots - hrom_snapshots)*100/ np.linalg.norm(rom_snapshots)
                    for i in range(len(mu_train)):
                        error = np.linalg.norm(rom_snapshots[:,i] - hrom_snapshots[:,i])*100/ np.linalg.norm(rom_snapshots[:,i])
                        mu_train_errors[i][3] = error
                    #######################################

            ##########################
            ###  Petrov Galerkin   ###
            elif chosen_projection_strategy == "petrov_galerkin":
                if any(item == "ROM" for item in training_stages):
                    fom_snapshots, n_clusters = self.__LaunchTrainROM(mu_train, use_non_linear_steps, cluster_list, clustering_method, n_clusters, bisection_tolerance,  POD_tolerance, use_aux_step)
                    # self._ChangeRomFlags(simulation_to_run = "TrainPG")
                    self.__LaunchTrainPG(mu_train, n_clusters, cluster_list, "TrainPG")
                    # self._ChangeRomFlags(simulation_to_run = "PG")
                    rom_snapshots = self.__LaunchROM(mu_train, cluster_list, "PG")
                    self.ROMvsFOM_train = np.linalg.norm(fom_snapshots - rom_snapshots)*100/ np.linalg.norm(fom_snapshots)
                    for i in range(len(mu_train)):
                        error = np.linalg.norm(fom_snapshots[:,i] - rom_snapshots[:,i])*100/ np.linalg.norm(fom_snapshots[:,i])
                        mu_train_errors.append([mu_train[i][0],mu_train[i][1],error,-1])
                if any(item == "HROM" for item in training_stages):
                    #FIXME there will be an error if we only train HROM, but not ROM
                    # self._ChangeRomFlags(simulation_to_run = "trainHROMPetrovGalerkin")
                    self.__LaunchTrainHROM(mu_train, n_clusters, cluster_list, "trainHROMPetrovGalerkin")
                    # self._ChangeRomFlags(simulation_to_run = "runHROMPetrovGalerkin")
                    hrom_snapshots = self.__LaunchHROM(mu_train, cluster_list, "runHROMPetrovGalerkin")
                    self.ROMvsHROM_train = np.linalg.norm(rom_snapshots - hrom_snapshots)*100/ np.linalg.norm(rom_snapshots)
                    for i in range(len(mu_train)):
                        error = np.linalg.norm(rom_snapshots[:,i] - hrom_snapshots[:,i])*100/ np.linalg.norm(rom_snapshots[:,i])
                        mu_train_errors[i][3] = error
            ##########################
            else:
                err_msg = f'Provided projection strategy {chosen_projection_strategy} is not supported. Available options are \'galerkin\', \'lspg\' and \'petrov_galerkin\'.'
                raise Exception(err_msg)
            np.save("centroids_list", self.centroids_to_plot)
        return mu_train_errors


    def Test(self, mu_test=[None]):
        mu_test_errors = []
        if (len(mu_test) > 0):
            chosen_projection_strategy = self.general_rom_manager_parameters["projection_strategy"].GetString()
            testing_stages = self.general_rom_manager_parameters["rom_stages_to_test"].GetStringArray()
            #######################
            ######  Galerkin ######
            if chosen_projection_strategy == "galerkin":
                if any(item == "ROM" for item in testing_stages):
                    fom_snapshots = self.__LaunchTestFOM(mu_test)
                    # self._ChangeRomFlags(simulation_to_run = "GalerkinROM")
                    rom_snapshots = self.__LaunchTestROM(mu_test, "GalerkinROM")
                    self.ROMvsFOM_test = np.linalg.norm(fom_snapshots - rom_snapshots)*100/ np.linalg.norm(fom_snapshots)
                    for i in range(len(mu_test)):
                        error = np.linalg.norm(fom_snapshots[:,i] - rom_snapshots[:,i])*100/ np.linalg.norm(fom_snapshots[:,i])
                        mu_test_errors.append([mu_test[i][0],mu_test[i][1],error,-1])

                if any(item == "HROM" for item in testing_stages):
                    #FIXME there will be an error if we only test HROM, but not ROM
                    self._ChangeRomFlags(simulation_to_run = "runHROMGalerkin")
                    hrom_snapshots = self.__LaunchTestHROM(mu_test)
                    self.ROMvsHROM_test = np.linalg.norm(rom_snapshots - hrom_snapshots)*100/ np.linalg.norm(rom_snapshots)
                    for i in range(len(mu_test)):
                        error = np.linalg.norm(rom_snapshots[:,i] - hrom_snapshots[:,i])*100/ np.linalg.norm(rom_snapshots[:,i])
                        mu_test_errors[i][3] = error


            #######################################
            ##  Least-Squares Petrov Galerkin   ###
            elif chosen_projection_strategy == "lspg":
                if any(item == "ROM" for item in testing_stages):
                    fom_snapshots = self.__LaunchTestFOM(mu_test)
                    # self._ChangeRomFlags(simulation_to_run = "lspg")
                    rom_snapshots = self.__LaunchTestROM(mu_test, "lspg")
                    self.ROMvsFOM_test = np.linalg.norm(fom_snapshots - rom_snapshots)*100/ np.linalg.norm(fom_snapshots)
                    for i in range(len(mu_test)):
                        error = np.linalg.norm(fom_snapshots[:,i] - rom_snapshots[:,i])*100/ np.linalg.norm(fom_snapshots[:,i])
                        mu_test_errors.append([mu_test[i][0],mu_test[i][1],error,-1])
                if any(item == "HROM" for item in testing_stages):
                    self._ChangeRomFlags(simulation_to_run = "runHROMLSPG")
                    hrom_snapshots = self.__LaunchTestHROM(mu_test)
                    self.ROMvsHROM_test = np.linalg.norm(rom_snapshots - hrom_snapshots)*100/ np.linalg.norm(rom_snapshots)
                    for i in range(len(mu_test)):
                        error = np.linalg.norm(rom_snapshots[:,i] - hrom_snapshots[:,i])*100/ np.linalg.norm(rom_snapshots[:,i])
                        mu_test_errors[i][3] = error
            #######################################


            ##########################
            ###  Petrov Galerkin   ###
            elif chosen_projection_strategy == "petrov_galerkin":
                if any(item == "ROM" for item in testing_stages):
                    fom_snapshots = self.__LaunchTestFOM(mu_test)
                    self._ChangeRomFlags(simulation_to_run = "PG")
                    rom_snapshots = self.__LaunchTestROM(mu_test)
                    self.ROMvsFOM_test = np.linalg.norm(fom_snapshots - rom_snapshots)*100/ np.linalg.norm(fom_snapshots)
                    for i in range(len(mu_test)):
                        error = np.linalg.norm(fom_snapshots[:,i] - rom_snapshots[:,i])*100/ np.linalg.norm(fom_snapshots[:,i])
                        mu_test_errors.append([mu_test[i][0],mu_test[i][1],error,-1])
                if any(item == "HROM" for item in testing_stages):
                    #FIXME there will be an error if we only train HROM, but not ROM
                    self._ChangeRomFlags(simulation_to_run = "runHROMPetrovGalerkin")
                    hrom_snapshots = self.__LaunchTestHROM(mu_test)
                    self.ROMvsHROM_test = np.linalg.norm(rom_snapshots - hrom_snapshots)*100/ np.linalg.norm(rom_snapshots)
                    for i in range(len(mu_test)):
                        error = np.linalg.norm(rom_snapshots[:,i] - hrom_snapshots[:,i])*100/ np.linalg.norm(rom_snapshots[:,i])
                        mu_test_errors[i][3] = error
            ##########################
            else:
                err_msg = f'Provided projection strategy {chosen_projection_strategy} is not supported. Available options are \'galerkin\', \'lspg\' and \'petrov_galerkin\'.'
                raise Exception(err_msg)
        return mu_test_errors



    def RunFOM(self, mu_run=[None]):
        self.__LaunchRunFOM(mu_run)

    def RunROM(self, mu_run=[None]):
        chosen_projection_strategy = self.general_rom_manager_parameters["projection_strategy"].GetString()
        #######################
        ######  Galerkin ######
        if chosen_projection_strategy == "galerkin":
            self._ChangeRomFlags(simulation_to_run = "GalerkinROM")
        #######################################
        ##  Least-Squares Petrov Galerkin   ###
        elif chosen_projection_strategy == "lspg":
            self._ChangeRomFlags(simulation_to_run = "lspg")
        ##########################
        ###  Petrov Galerkin   ###
        elif chosen_projection_strategy == "petrov_galerkin":
            self._ChangeRomFlags(simulation_to_run = "PG")
        else:
            err_msg = f'Provided projection strategy {chosen_projection_strategy} is not supported. Available options are \'galerkin\', \'lspg\' and \'petrov_galerkin\'.'
            raise Exception(err_msg)
        self.__LaunchRunROM(mu_run)

    def RunHROM(self, mu_run=[None], use_full_model_part = False):
        chosen_projection_strategy = self.general_rom_manager_parameters["projection_strategy"].GetString()
        #######################
        ######  Galerkin ######
        if chosen_projection_strategy == "galerkin":
            self._ChangeRomFlags(simulation_to_run = "runHROMGalerkin")
        #######################################
        ##  Least-Squares Petrov Galerkin   ###
        elif chosen_projection_strategy == "lspg":
            self._ChangeRomFlags(simulation_to_run = "runHROMLSPG")
        ##########################
        ###  Petrov Galerkin   ###
        elif chosen_projection_strategy == "petrov_galerkin":
            self._ChangeRomFlags(simulation_to_run = "runHROMPetrovGalerkin")
        else:
            err_msg = f'Provided projection strategy {chosen_projection_strategy} is not supported. Available options are \'galerkin\', \'lspg\' and \'petrov_galerkin\'.'
            raise Exception(err_msg)
        self.__LaunchRunHROM(mu_run, use_full_model_part)


    def PrintErrors(self):
        training_stages = self.general_rom_manager_parameters["rom_stages_to_train"].GetStringArray()
        testing_stages = self.general_rom_manager_parameters["rom_stages_to_test"].GetStringArray()
        if any(item == "ROM" for item in training_stages) and self.ROMvsFOM_train != 0:
            print("approximation error in train set FOM vs ROM: ", self.ROMvsFOM_train," %")
        if any(item == "HROM" for item in training_stages) and self.ROMvsHROM_train != 0:
            print("approximation error in train set ROM vs HROM: ", self.ROMvsHROM_train," %")
            print("average approximation error per node in train set ROM vs HROM: ", self.ROMvsHROM_train," %")
        if any(item == "ROM" for item in testing_stages) and self.ROMvsFOM_test != 0:
            print("approximation error in test set FOM vs ROM: ", self.ROMvsFOM_test," %")
        if any(item == "HROM" for item in testing_stages) and self.ROMvsHROM_test != 0:
            print("approximation error in test set ROM vs HROM: ",  self.ROMvsHROM_test," %")
            print("average approximation error per node in test set ROM vs HROM: ",  self.ROMvsHROM_test," %")

    def __LaunchTrainROM(self, mu_train, use_non_linear_steps, cluster_list, clustering_method, n_clusters, bisection_tolerance,  POD_tolerance, use_aux_step):
        """
        This method should be parallel capable
        """
        if not os.path.exists("DataBase"):
            os.mkdir("DataBase")
            os.mkdir("DataBase/Results")
            os.mkdir("DataBase/Data")
            os.mkdir("DataBase/NonLinearStepsData")
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        # Build data base and compute the clusters:
        aux_SnapshotsMatrix = []; full_SnapshotsMatrix = []
        for Id, mu in enumerate(mu_train):
            parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
            parameters_copy = self._AddBasisCreationToProjectParameters(parameters_copy) #TODO stop using the RomBasisOutputProcess to store the snapshots. Use instead the upcoming build-in function
            parameters_copy = self._StoreResultsByName(parameters_copy,'FOM_Fit',mu,Id)
            parameters_copy["output_processes"]["vtk_output"][0]["Parameters"]["output_path"].SetString("DataBase/Results/"+ str(mu[0]) + ", " + str(mu[1]))
            model = KratosMultiphysics.Model()
            analysis_stage_class = self._GetAnalysisStageClass(parameters_copy)
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy)
            case_name = "DataBase/"+ str(mu[0]) + ", " + str(mu[1]) +".npy"
            non_linear_steps_case_name = "DataBase/NonLinearStepsData/"+ str(mu[0]) + ", " + str(mu[1]) +".npy"
            if not os.path.exists(case_name) or not os.path.exists(non_linear_steps_case_name):
                simulation.Run()
                for process in simulation._GetListOfOutputProcesses():
                    if isinstance(process, CalculateRomBasisOutputProcess):
                        BasisOutputProcess = process
                np.save(non_linear_steps_case_name, np.array(simulation._GetSolver()._GetSolutionStrategy().GetIntermediateSolutionsMatrix()))
                np.save(case_name, BasisOutputProcess._GetSnapshotsMatrix())
                aux_SnapshotsMatrix.append(BasisOutputProcess._GetSnapshotsMatrix()) #TODO add a CustomMethod() as a standard method in the Analysis Stage to retrive some solution
                full_SnapshotsMatrix.append(BasisOutputProcess._GetSnapshotsMatrix())
                full_SnapshotsMatrix.append(np.array(simulation._GetSolver()._GetSolutionStrategy().GetIntermediateSolutionsMatrix()))
            else:
                case = np.load(case_name)
                non_linear_steps_case = np.load(non_linear_steps_case_name)
                simulation.Initialize()
                for process in simulation._GetListOfOutputProcesses():
                    if isinstance(process, CalculateRomBasisOutputProcess):
                        BasisOutputProcess = process
                aux_SnapshotsMatrix.append(case)
                full_SnapshotsMatrix.append(case)
                full_SnapshotsMatrix.append(non_linear_steps_case)
        aux_SnapshotsMatrix = np.block(aux_SnapshotsMatrix)
        full_SnapshotsMatrix = np.block(full_SnapshotsMatrix)

        #launch clustering
        if clustering_method=="k-means":
             n_clusters = self.kmeans_method(aux_SnapshotsMatrix.T, n_clusters, np.array(np.double(mu_train)), cluster_list)
        elif clustering_method=="PEBL":
            n_clusters = self.pebl_method(aux_SnapshotsMatrix.T, bisection_tolerance,  POD_tolerance, cluster_list, list=mu_train)
        elif clustering_method=="PEBL_full":
            n_clusters = self.pebl_method_full_snapshots_matrix(full_SnapshotsMatrix. T,aux_SnapshotsMatrix.T, bisection_tolerance,  POD_tolerance, cluster_list, list=mu_train)
        else:
            raise RuntimeError("'" + clustering_method + "' method doesn't exist.")
        
        NumberOfRomModes_list = np.zeros(n_clusters)
        for i in range(n_clusters):
            SnapshotsMatrix = []
            mu_aux = [mu for mu, cluster in zip(mu_train, cluster_list[:][0]) if cluster == i]
            for Id, mu in enumerate(mu_aux):
                parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
                parameters_copy = self._AddBasisCreationToProjectParameters(parameters_copy) #TODO stop using the RomBasisOutputProcess to store the snapshots. Use instead the upcoming build-in function
                parameters_copy = self._StoreResultsByName(parameters_copy,'FOM_Fit',mu,Id)
                parameters_copy["output_processes"]["rom_output"][0]["Parameters"]["rom_basis_output_folder"].SetString("rom_data_cluster_"+str(i))
                parameters_copy["output_processes"]["vtk_output"][0]["Parameters"]["output_path"].SetString("DataBase/Results/"+ str(mu[0]) + ", " + str(mu[1]))
                model = KratosMultiphysics.Model()
                analysis_stage_class = self._GetAnalysisStageClass(parameters_copy)
                simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy)
                if use_non_linear_steps:
                    case_name = "DataBase/"+ str(mu[0]) + ", " + str(mu[1]) +".npy"
                    non_linear_steps_case_name = "DataBase/NonLinearStepsData/"+ str(mu[0]) + ", " + str(mu[1]) +".npy"
                    case = np.load(case_name)
                    non_linear_steps_case = np.load(non_linear_steps_case_name)
                    simulation.Initialize()
                    for process in simulation._GetListOfOutputProcesses():
                        if isinstance(process, CalculateRomBasisOutputProcess):
                            BasisOutputProcess = process
                    SnapshotsMatrix.append(case)
                    SnapshotsMatrix.append(non_linear_steps_case) 
                else:
                    case_name = "DataBase/"+ str(mu[0]) + ", " + str(mu[1]) +".npy"
                    case = np.load(case_name)
                    simulation.Initialize()
                    for process in simulation._GetListOfOutputProcesses():
                        if isinstance(process, CalculateRomBasisOutputProcess):
                            BasisOutputProcess = process
                    SnapshotsMatrix.append(case)
                if use_aux_step:
                    aux_step = np.load("DataBase/aux_data/"+ str(mu[0]) + ", " + str(mu[1]) +".npy")
                    if aux_step.ndim ==1: #check if matrix contains a single mode (a 1D numpy array)
                        aux_step.reshape(-1,1)
                    SnapshotsMatrix.append(aux_step) 
            if len(mu_aux)>0:
                SnapshotsMatrix = np.block(SnapshotsMatrix)
                BasisOutputProcess._PrintRomBasis(SnapshotsMatrix) #Calling the RomOutput Process for creating the RomParameter.json
                NumberOfRomModes_list[i] = BasisOutputProcess.GetNumberOfRomModes()
        self.PlotClusterListWithRomModes(n_clusters,mu_train,cluster_list,NumberOfRomModes_list,clustering_method)
        return aux_SnapshotsMatrix, n_clusters


    def __LaunchROM(self, mu_train, cluster_list, simulation_to_run):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        SnapshotsMatrix = []
        for Id, mu in enumerate(mu_train):
            self.general_rom_manager_parameters["ROM"]["rom_basis_output_folder"].SetString("rom_data_cluster_"+str(np.int0(cluster_list[0][Id])))
            self._ChangeRomFlags(simulation_to_run)
            parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
            parameters_copy = self._AddBasisCreationToProjectParameters(parameters_copy)  #TODO stop using the RomBasisOutputProcess to store the snapshots. Use instead the upcoming build-in function
            parameters_copy = self._StoreResultsByName(parameters_copy,'ROM_Fit',mu,Id)
            # parameters_copy["output_processes"]["rom_output"][0]["Parameters"]["rom_basis_output_folder"].SetString("rom_data_cluster_"+str(np.int0(cluster_list[0][Id])))
            model = KratosMultiphysics.Model()
            analysis_stage_class = type(SetUpSimulationInstance(model, parameters_copy))
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy)
            simulation.Run()
            model_part = simulation.model[parameters_copy["solver_settings"]["model_part_name"].GetString()].GetRootModelPart()
            self.train_info_steps_list.append([Id, 
                                         simulation._GetSolver()._GetBuilderAndSolver().GetCorrectionsCounter(),
                                         model_part.ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER],
                                         model_part.ProcessInfo[KratosMultiphysics.CONVERGENCE_RATIO],
                                         ])
            for process in simulation._GetListOfOutputProcesses():
                if isinstance(process, CalculateRomBasisOutputProcess):
                    BasisOutputProcess = process
            SnapshotsMatrix.append(BasisOutputProcess._GetSnapshotsMatrix()) #TODO add a CustomMethod() as a standard method in the Analysis Stage to retrive some solution
        SnapshotsMatrix = np.block(SnapshotsMatrix)
        print("Training simulations info:")
        print("Id, N Corrections, NL iterations, Convergence ratio")
        print(self.train_info_steps_list)

        import openpyxl
        wb = openpyxl.Workbook()
        hoja = wb.active
        hoja.append(('Id', 'N Corrections', 'NL iterations', 'Convergence ratio'))
        for item in self.train_info_steps_list:
            hoja.append(item)
        wb.save('train_data.xlsx')

        return SnapshotsMatrix


    def __LaunchTrainPG(self, mu_train, n_clusters, cluster_list, simulation_to_run):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        for i in range(n_clusters):
            PetrovGalerkinTrainMatrix = []
            mu_aux = [mu for mu, cluster in zip(mu_train, cluster_list[:][0]) if cluster == i]
            for Id, mu in enumerate(mu_aux):
                self.general_rom_manager_parameters["ROM"]["rom_basis_output_folder"].SetString("rom_data_cluster_"+str(i))
                self._ChangeRomFlags(simulation_to_run)
                parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
                parameters_copy = self._AddBasisCreationToProjectParameters(parameters_copy) 
                parameters_copy = self._StoreNoResults(parameters_copy)
                # parameters_copy["output_processes"]["rom_output"][0]["Parameters"]["rom_basis_output_folder"].SetString("rom_data_cluster_"+str(i))
                model = KratosMultiphysics.Model()
                analysis_stage_class = type(SetUpSimulationInstance(model, parameters_copy))
                simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy)
                simulation.Run()
                PetrovGalerkinTrainMatrix.append(simulation.GetPetrovGalerkinTrainUtility()._GetSnapshotsMatrix()) #TODO is the best way of extracting the Projected Residuals calling the HROM residuals utility?
                # PetrovGalerkinTrainMatrix.append(np.array(simulation._GetSolver()._GetSolutionStrategy().GetIntermediateSolutionsMatrix()))
            simulation.GetPetrovGalerkinTrainUtility().CalculateAndSaveBasis(np.block(PetrovGalerkinTrainMatrix))


    def __LaunchTrainHROM(self, mu_train, n_clusters, cluster_list, simulation_to_run):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        for i in range(n_clusters):
            RedidualsSnapshotsMatrix = []
            mu_aux = [mu for mu, cluster in zip(mu_train, cluster_list[:][0]) if cluster == i]
            for Id, mu in enumerate(mu_aux):
                self.general_rom_manager_parameters["ROM"]["rom_basis_output_folder"].SetString("rom_data_cluster_"+str(i))
                self._ChangeRomFlags(simulation_to_run)
                parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
                parameters_copy = self._AddBasisCreationToProjectParameters(parameters_copy)
                parameters_copy = self._StoreNoResults(parameters_copy)
                # parameters_copy["output_processes"]["rom_output"][0]["Parameters"]["rom_basis_output_folder"].SetString("rom_data_cluster_"+str(i))
                model = KratosMultiphysics.Model()
                analysis_stage_class = type(SetUpSimulationInstance(model, parameters_copy))
                simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy)
                simulation.Run()
                RedidualsSnapshotsMatrix.append(simulation.GetHROM_utility()._GetResidualsProjectedMatrix()) #TODO is the best way of extracting the Projected Residuals calling the HROM residuals utility?
            RedidualsSnapshotsMatrix = np.block(RedidualsSnapshotsMatrix)
            u,_,_,_ = RandomizedSingularValueDecomposition(COMPUTE_V=False).Calculate(RedidualsSnapshotsMatrix,
            self.SetHromTrainingParameters()["element_selection_svd_truncation_tolerance"].GetDouble())
            simulation.GetHROM_utility().hyper_reduction_element_selector.SetUp(u, InitialCandidatesSet = simulation.GetHROM_utility().candidate_ids)
            simulation.GetHROM_utility().hyper_reduction_element_selector.Run()
            if not simulation.GetHROM_utility().hyper_reduction_element_selector.success:
                KratosMultiphysics.Logger.PrintWarning("HRomTrainingUtility", "The Empirical Cubature Method did not converge using the initial set of candidates. Launching again without initial candidates.")
                #Imposing an initial candidate set can lead to no convergence. Restart without imposing the initial candidate set
                simulation.GetHROM_utility().hyper_reduction_element_selector.SetUp(u, InitialCandidatesSet = None)
                simulation.GetHROM_utility().hyper_reduction_element_selector.Run()
            simulation.GetHROM_utility().AppendHRomWeightsToRomParameters()
            # simulation.GetHROM_utility().CreateHRomModelParts()


    def __LaunchHROM(self, mu_train, cluster_list, simulation_to_run):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        SnapshotsMatrix = []
        for Id, mu in enumerate(mu_train):
            self.general_rom_manager_parameters["ROM"]["rom_basis_output_folder"].SetString("rom_data_cluster_"+str(np.int0(cluster_list[0][Id])))
            self._ChangeRomFlags(simulation_to_run)
            parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
            parameters_copy = self._AddBasisCreationToProjectParameters(parameters_copy)
            parameters_copy = self._StoreResultsByName(parameters_copy,'HROM_Fit',mu,Id)
            # parameters_copy["output_processes"]["rom_output"][0]["Parameters"]["rom_basis_output_folder"].SetString("rom_data_cluster_"+str(np.int0(cluster_list[0][Id])))
            model = KratosMultiphysics.Model()
            analysis_stage_class = type(SetUpSimulationInstance(model, parameters_copy))
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy)
            simulation.Run()
            for process in simulation._GetListOfOutputProcesses():
                if isinstance(process, CalculateRomBasisOutputProcess):
                    BasisOutputProcess = process
            SnapshotsMatrix.append(BasisOutputProcess._GetSnapshotsMatrix()) #TODO add a CustomMethod() as a standard method in the Analysis Stage to retrive some solution
        SnapshotsMatrix = np.block(SnapshotsMatrix)

        return SnapshotsMatrix


    def __LaunchTestFOM(self, mu_test):
        """
        This method should be parallel capable
        """
        # FIXME We must include a method to retrive solutions in the nodes and stop using the CalculateRomBasisOutputProcess to stote snapshots
        if not os.path.exists("DataBase"):
            os.mkdir("DataBase")
            os.mkdir("DataBase/Results")
            os.mkdir("DataBase/Data")
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        SnapshotsMatrix = []
        for Id, mu in enumerate(mu_test):
            parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
            parameters_copy = self._AddBasisCreationToProjectParameters(parameters_copy) #TODO stop using the RomBasisOutputProcess to store the snapshots. Use instead the upcoming build-in function
            parameters_copy = self._StoreResultsByName(parameters_copy,'FOM_Test',mu,Id)
            model = KratosMultiphysics.Model()
            analysis_stage_class = self._GetAnalysisStageClass(parameters_copy)
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy)
            case_name = "DataBase/"+ str(mu[0]) + ", " + str(mu[1]) +".npy"
            if not os.path.exists(case_name):
                simulation.Run()
                for process in simulation._GetListOfOutputProcesses():
                    if isinstance(process, CalculateRomBasisOutputProcess):
                        BasisOutputProcess = process
                np.save(case_name, BasisOutputProcess._GetSnapshotsMatrix())
                SnapshotsMatrix.append(BasisOutputProcess._GetSnapshotsMatrix()) #TODO add a CustomMethod() as a standard method in the Analysis Stage to retrive some solution
            else:
                case = np.load(case_name)
                simulation.Initialize()
                for process in simulation._GetListOfOutputProcesses():
                    if isinstance(process, CalculateRomBasisOutputProcess):
                        BasisOutputProcess = process
                SnapshotsMatrix.append(case)
        SnapshotsMatrix = np.block(SnapshotsMatrix)
        return SnapshotsMatrix


    def __LaunchTestROM(self, mu_test, simulation_to_run):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        SnapshotsMatrix = []
        cluster_list = np.zeros(len(mu_test))
        centroids_list = np.load("centroids_list.npy")        
        for Id, mu in enumerate(mu_test):
            cluster_list[Id] = self.GetClusterId(centroids_list, [mu[1],mu[0]])
            self.general_rom_manager_parameters["ROM"]["rom_basis_output_folder"].SetString("rom_data_cluster_"+str(self.GetClusterId(centroids_list, [mu[1],mu[0]])))
            self._ChangeRomFlags(simulation_to_run)
            parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
            parameters_copy = self._AddBasisCreationToProjectParameters(parameters_copy)  #TODO stop using the RomBasisOutputProcess to store the snapshots. Use instead the upcoming build-in function
            parameters_copy = self._StoreResultsByName(parameters_copy,'ROM_Test',mu,Id)
            model = KratosMultiphysics.Model()
            analysis_stage_class = type(SetUpSimulationInstance(model, parameters_copy))
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy)
            simulation.Run()
            model_part = simulation.model[parameters_copy["solver_settings"]["model_part_name"].GetString()].GetRootModelPart()
            self.test_info_steps_list.append([Id, 
                                         simulation._GetSolver()._GetBuilderAndSolver().GetCorrectionsCounter(),
                                         model_part.ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER],
                                         model_part.ProcessInfo[KratosMultiphysics.CONVERGENCE_RATIO],
                                         ])         
            for process in simulation._GetListOfOutputProcesses():
                if isinstance(process, CalculateRomBasisOutputProcess):
                    BasisOutputProcess = process
            SnapshotsMatrix.append(BasisOutputProcess._GetSnapshotsMatrix()) #TODO add a CustomMethod() as a standard method in the Analysis Stage to retrive some solution
        SnapshotsMatrix = np.block(SnapshotsMatrix)

        print("Training simulations info:")
        print("Id, N Corrections, NL iterations, Convergence ratio")
        print(self.train_info_steps_list)

        print("Testing simulations info:")
        print("Id, N Corrections, NL iterations, Convergence ratio")
        print(self.test_info_steps_list)

        #plot test values
        fig = plt.figure()
        fig.set_figwidth(15.0)
        fig.set_figheight(10.0)
        for i in range(len(centroids_list)):
            mu_aux = [mu for mu, cluster in zip(mu_test, cluster_list) if cluster == i]
            mu_train_a = np.zeros(len(mu_aux))
            mu_train_m = np.zeros(len(mu_aux))
            for j in range(len(mu_aux)):
                mu_train_a[j] = mu_aux[j][0]
                mu_train_m[j] = mu_aux[j][1]
            fig = plt.scatter(mu_train_m, mu_train_a, label="Cluster " + str(i))
            fig = plt.scatter(centroids_list[i,0], centroids_list[i,1], c='k',marker="s", s= 150)
            fig = plt.text(centroids_list[i,0], centroids_list[i,1]-0.05, "Centroid C"+str(i), ha='center')
            for j in range(len(mu_aux)):
                fig = plt.text(mu_train_m[j], mu_train_a[j]-0.05, str(i), ha='center')
        fig = plt.title(" Train values and clusters")
        fig = plt.ylabel('Alpha')
        fig = plt.xlabel('Mach')
        fig = plt.legend(bbox_to_anchor=(1.004, 0.9, 1.0, 0.102), loc='upper left', borderaxespad=0.)
        fig = plt.savefig('mu test and clusters.png',bbox_inches='tight' )
        fig = plt.close('all')

        import openpyxl
        wb = openpyxl.Workbook()
        hoja = wb.active
        hoja.append(('Id', 'N Corrections', 'NL iterations', 'Convergence ratio'))
        for item in self.test_info_steps_list:
            hoja.append(item)
        wb.save('test_data.xlsx')

        return SnapshotsMatrix



    def __LaunchTestHROM(self, mu_test):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        SnapshotsMatrix = []
        for Id, mu in enumerate(mu_test):
            parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
            parameters_copy = self._AddBasisCreationToProjectParameters(parameters_copy)
            parameters_copy = self._StoreResultsByName(parameters_copy,'HROM_Test',mu,Id)
            model = KratosMultiphysics.Model()
            analysis_stage_class = type(SetUpSimulationInstance(model, parameters_copy))
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy)
            simulation.Run()
            for process in simulation._GetListOfOutputProcesses():
                if isinstance(process, CalculateRomBasisOutputProcess):
                    BasisOutputProcess = process
            SnapshotsMatrix.append(BasisOutputProcess._GetSnapshotsMatrix()) #TODO add a CustomMethod() as a standard method in the Analysis Stage to retrive some solution
        SnapshotsMatrix = np.block(SnapshotsMatrix)

        return SnapshotsMatrix


    def __LaunchRunFOM(self, mu_run):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        for Id, mu in enumerate(mu_run):
            parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
            parameters_copy = self._StoreResultsByName(parameters_copy,'FOM_Run',mu,Id)
            model = KratosMultiphysics.Model()
            analysis_stage_class = self._GetAnalysisStageClass(parameters_copy)
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy)
            simulation.Run()


    def __LaunchRunROM(self, mu_run):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        for Id, mu in enumerate(mu_run):
            parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
            parameters_copy = self._StoreResultsByName(parameters_copy,'ROM_Run',mu,Id)
            model = KratosMultiphysics.Model()
            analysis_stage_class = type(SetUpSimulationInstance(model, parameters_copy))
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy)
            simulation.Run()


    def __LaunchRunHROM(self, mu_run, use_full_model_part):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        if not use_full_model_part:
            model_part_name = parameters["solver_settings"]["model_import_settings"]["input_filename"].GetString()
            parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(f"{model_part_name}HROM")

        for Id, mu in enumerate(mu_run):
            parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
            parameters_copy = self._StoreResultsByName(parameters_copy,'HROM_Run',mu,Id)
            model = KratosMultiphysics.Model()
            analysis_stage_class = type(SetUpSimulationInstance(model, parameters_copy))
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy)
            simulation.Run()


    def _AddHromParametersToRomParameters(self,f):
        f["hrom_settings"]["element_selection_type"] = self.SetHromTrainingParameters()["element_selection_type"].GetString()
        f["hrom_settings"]["element_selection_svd_truncation_tolerance"] = self.SetHromTrainingParameters()["element_selection_svd_truncation_tolerance"].GetDouble()
        f["hrom_settings"]["create_hrom_visualization_model_part"] = self.SetHromTrainingParameters()["create_hrom_visualization_model_part"].GetBool()
        f["hrom_settings"]["echo_level"] = self.SetHromTrainingParameters()["echo_level"].GetInt()
        f["hrom_settings"]["include_condition_parents"] = self.SetHromTrainingParameters()["include_condition_parents"].GetBool()
        f["hrom_settings"]["initial_candidate_elements_model_part_list"] = self.SetHromTrainingParameters()["initial_candidate_elements_model_part_list"].GetStringArray()
        f["hrom_settings"]["initial_candidate_conditions_model_part_list"] = self.SetHromTrainingParameters()["initial_candidate_conditions_model_part_list"].GetStringArray()
        f["hrom_settings"]["constraint_sum_weights"] = self.SetHromTrainingParameters()["constraint_sum_weights"].GetBool()

    def _ChangeRomFlags(self, simulation_to_run = 'ROM'):
        """
        This method updates the Flags present in the RomParameters.json file
        for launching the correct part of the ROM workflow
        """
        #other options: "trainHROM", "runHROM"
        parameters_file_folder = self.general_rom_manager_parameters["ROM"]["rom_basis_output_folder"].GetString() if self.general_rom_manager_parameters["ROM"].Has("rom_basis_output_folder") else "rom_data"
        parameters_file_name = self.general_rom_manager_parameters["ROM"]["rom_basis_output_name"].GetString() if self.general_rom_manager_parameters["ROM"].Has("rom_basis_output_name") else "RomParameters"

        # Convert to Path objects
        parameters_file_folder = Path(parameters_file_folder)
        parameters_file_name = Path(parameters_file_name)

        parameters_file_path = parameters_file_folder / parameters_file_name.with_suffix('.json')

        with parameters_file_path.open('r+') as parameter_file:
            f=json.load(parameter_file)
            f['assembling_strategy'] = self.general_rom_manager_parameters['assembling_strategy'].GetString() if self.general_rom_manager_parameters.Has('assembling_strategy') else 'global'
            self._AddHromParametersToRomParameters(f)
            if simulation_to_run=='GalerkinROM':
                f['projection_strategy']="galerkin"
                f['train_hrom']=False
                f['run_hrom']=False
                f["rom_settings"]['rom_bns_settings'] = self._SetGalerkinBnSParameters()
            elif simulation_to_run=='trainHROMGalerkin':
                f['train_hrom']=True
                f['run_hrom']=False
                f["rom_settings"]['rom_bns_settings'] = self._SetGalerkinBnSParameters()
            elif simulation_to_run=='runHROMGalerkin':
                f['projection_strategy']="galerkin"
                f['train_hrom']=False
                f['run_hrom']=True
                f["rom_settings"]['rom_bns_settings'] = self._SetGalerkinBnSParameters()
            elif simulation_to_run == 'lspg':
                f['train_hrom'] = False
                f['run_hrom'] = False
                f['projection_strategy'] = "lspg"
                f["rom_settings"]['rom_bns_settings'] = self._SetLSPGBnSParameters()
            elif simulation_to_run == 'trainHROMLSPG':
                f['train_hrom'] = True
                f['run_hrom'] = False
                f['projection_strategy'] = "lspg"
                f["rom_settings"]['rom_bns_settings'] = self._SetLSPGBnSParameters()
            elif simulation_to_run == 'runHROMLSPG':
                f['train_hrom'] = False
                f['run_hrom'] = True
                f['projection_strategy'] = "lspg"
                f["rom_settings"]['rom_bns_settings'] = self._SetLSPGBnSParameters()
            elif simulation_to_run == 'TrainPG':
                f['train_hrom'] = False
                f['run_hrom'] = False
                f['projection_strategy'] = "lspg"
                f["rom_settings"]['rom_bns_settings'] = self._SetLSPGBnSParameters()
                f["rom_settings"]['rom_bns_settings']['train_petrov_galerkin'] = True  # Override the default
            elif simulation_to_run=='PG':
                f['train_hrom']=False
                f['run_hrom']=False
                f['projection_strategy']="petrov_galerkin"
                f["rom_settings"]['rom_bns_settings'] = self._SetPetrovGalerkinBnSParameters()
            elif simulation_to_run=='trainHROMPetrovGalerkin':
                f['train_hrom']=True
                f['run_hrom']=False
                f['projection_strategy']="petrov_galerkin"
                f["rom_settings"]['rom_bns_settings'] = self._SetPetrovGalerkinBnSParameters()
            elif simulation_to_run=='runHROMPetrovGalerkin':
                f['train_hrom']=False
                f['run_hrom']=True
                f['projection_strategy']="petrov_galerkin"
                f["rom_settings"]['rom_bns_settings'] = self._SetPetrovGalerkinBnSParameters()
            else:
                raise Exception(f'Unknown flag "{simulation_to_run}" change for RomParameters.json')
            parameter_file.seek(0)
            json.dump(f,parameter_file,indent=4)
            parameter_file.truncate()

    def _SetGalerkinBnSParameters(self):
        # Retrieve the default parameters as a JSON string and parse it into a dictionary
        defaults_json = self._GetGalerkinBnSParameters()
        defaults = json.loads(defaults_json)

        # Ensure 'galerkin_rom_bns_settings' exists in ROM parameters
        if not self.general_rom_manager_parameters["ROM"].Has("galerkin_rom_bns_settings"):
            self.general_rom_manager_parameters["ROM"].AddEmptyValue("galerkin_rom_bns_settings")

        # Get the ROM parameters for Galerkin
        rom_params = self.general_rom_manager_parameters["ROM"]["galerkin_rom_bns_settings"]

        # Update defaults with any existing ROM parameters
        self._UpdateDefaultsWithRomParams(defaults, rom_params)

        return defaults

    def _SetLSPGBnSParameters(self):
        # Retrieve the default parameters as a JSON string and parse it into a dictionary
        defaults_json = self._GetLSPGBnSParameters()
        defaults = json.loads(defaults_json)

        # Ensure 'lspg_rom_bns_settings' exists in ROM parameters
        if not self.general_rom_manager_parameters["ROM"].Has("lspg_rom_bns_settings"):
            self.general_rom_manager_parameters["ROM"].AddEmptyValue("lspg_rom_bns_settings")

        # Get the ROM parameters for LSPG
        rom_params = self.general_rom_manager_parameters["ROM"]["lspg_rom_bns_settings"]

        # Update defaults with any existing ROM parameters
        self._UpdateDefaultsWithRomParams(defaults, rom_params)

        return defaults

    def _SetPetrovGalerkinBnSParameters(self):
        # Retrieve the default parameters as a JSON string and parse it into a dictionary
        defaults_json = self._GetPetrovGalerkinBnSParameters()
        defaults = json.loads(defaults_json)

        # Ensure 'petrov_galerkin_rom_bns_settings' exists in ROM parameters
        if not self.general_rom_manager_parameters["ROM"].Has("petrov_galerkin_rom_bns_settings"):
            self.general_rom_manager_parameters["ROM"].AddEmptyValue("petrov_galerkin_rom_bns_settings")

        # Get the ROM parameters for Petrov-Galerkin
        rom_params = self.general_rom_manager_parameters["ROM"]["petrov_galerkin_rom_bns_settings"]

        # Update defaults with any existing ROM parameters
        self._UpdateDefaultsWithRomParams(defaults, rom_params)

        return defaults

    def _UpdateDefaultsWithRomParams(self, defaults, rom_params):
        for key, default_value in defaults.items():
            if rom_params.Has(key):
                if isinstance(default_value, bool):
                    defaults[key] = rom_params[key].GetBool()
                elif isinstance(default_value, str):
                    defaults[key] = rom_params[key].GetString()
                elif isinstance(default_value, float):
                    defaults[key] = rom_params[key].GetDouble()
        return defaults

    def _AddBasisCreationToProjectParameters(self, parameters):
        #FIXME make sure no other rom_output already existed. If so, erase the prior and keep only the one in self._SetRomTrainingParameters()
        parameters["output_processes"].AddEmptyArray("rom_output")
        parameters["output_processes"]["rom_output"].Append(self._SetRomTrainingParameters())

        return parameters




    def _StoreResultsByName(self,parameters,results_name,mu, Id):

        if  self.general_rom_manager_parameters["output_name"].GetString() == "mu":
            case_name = (", ".join(map(str, mu)))
        elif self.general_rom_manager_parameters["output_name"].GetString() == "id":
            case_name = str(Id)
        if self.general_rom_manager_parameters["save_gid_output"].GetBool():
            parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString('Results/'+ results_name +  case_name)
        else:
            parameters["output_processes"].RemoveValue("gid_output")
        if self.general_rom_manager_parameters["save_vtk_output"].GetBool():
            parameters["output_processes"]["vtk_output"][0]["Parameters"]["output_path"].SetString('Results/'+ results_name + case_name)
        else:
            parameters["output_processes"].RemoveValue("vtk_output")

        return parameters


    def _StoreNoResults(self, parameters):
        parameters["output_processes"].RemoveValue("gid_output")
        parameters["output_processes"].RemoveValue("vtk_output")

        return parameters



    def _SetRomTrainingParameters(self):
        defaults = self._GetDefaulRomBasisOutputParameters()
        defaults["Parameters"]["rom_manager"].SetBool(True)  # Set the flag to true when inside the RomManager to trigger particular behavior for multiple parameters

        rom_params = self.general_rom_manager_parameters["ROM"]

        keys_to_copy = [
            "svd_truncation_tolerance",
            "model_part_name",
            "rom_basis_output_format",
            "rom_basis_output_name",
            "rom_basis_output_folder",
            "nodal_unknowns",
            "snapshots_interval",
        ]

        for key in keys_to_copy:
            if key in rom_params.keys():
                defaults["Parameters"][key] = rom_params[key]

        return defaults




    def SetHromTrainingParameters(self):
        defaults = self._GetDefaulHromTrainingParameters()
        hrom_params = self.general_rom_manager_parameters["HROM"]

        keys_to_copy = [
            "element_selection_type",
            "element_selection_svd_truncation_tolerance",
            "create_hrom_visualization_model_part",
            "echo_level",
            "constraint_sum_weights",
        ]

        for key in keys_to_copy:
            if key in hrom_params.keys():
                defaults[key] = hrom_params[key]

        return defaults



    def _GetAnalysisStageClass(self, parameters):

        analysis_stage_module_name = parameters["analysis_stage"].GetString()
        analysis_stage_class_name = analysis_stage_module_name.split('.')[-1]
        analysis_stage_class_name = ''.join(x.title() for x in analysis_stage_class_name.split('_'))

        analysis_stage_module = importlib.import_module(analysis_stage_module_name)
        analysis_stage_class = getattr(analysis_stage_module, analysis_stage_class_name)

        return analysis_stage_class


    def _GetDefaulRomBasisOutputParameters(self):
        rom_training_parameters = KratosMultiphysics.Parameters("""{
                "python_module" : "calculate_rom_basis_output_process",
                "kratos_module" : "KratosMultiphysics.RomApplication",
                "process_name"  : "CalculateRomBasisOutputProcess",
                "help"          : "This process should write the Rom basis",
                "Parameters"    :
                {
                    "model_part_name": "",
                    "rom_manager" : false,      // set to false for manual manipulation of ROM via flags in the RomParameters
                    "snapshots_control_type": "step",
                    "snapshots_interval": 1.0,
                    "nodal_unknowns":  [],
                    "rom_basis_output_format": "json",
                    "rom_basis_output_name": "RomParameters",
                    "rom_basis_output_folder": "rom_data",
                    "svd_truncation_tolerance": 1e-3
                }
            }""")
        return rom_training_parameters


    def _GetDefaulHromTrainingParameters(self):
        hrom_training_parameters = KratosMultiphysics.Parameters("""{
                "hrom_format": "numpy",
                "element_selection_type": "empirical_cubature",
                "element_selection_svd_truncation_tolerance": 1.0e-6,
                "echo_level" : 0,
                "create_hrom_visualization_model_part" : true,
                "projection_strategy": "galerkin",
                "include_conditions_model_parts_list": [],
                "include_elements_model_parts_list": [],
                "initial_candidate_elements_model_part_list" : [],
                "initial_candidate_conditions_model_part_list" : [],
                "include_nodal_neighbouring_elements_model_parts_list":[],
                "include_minimum_condition": false,
                "include_condition_parents": false,
                "constraint_sum_weights": true
            }""")
        return hrom_training_parameters


    def _GetGalerkinBnSParameters(self):
        # Define the default settings in JSON format for Galerkin BnS
        rom_bns_settings = """{
            "monotonicity_preserving": false
        }"""
        return rom_bns_settings

    def _GetPetrovGalerkinBnSParameters(self):
        # Define the default settings in JSON format for Petrov-Galerkin BnS
        rom_bns_settings = """{
            "monotonicity_preserving": false
        }"""
        return rom_bns_settings

    def _GetLSPGBnSParameters(self):
        # Define the default settings in JSON format
        # Comments:
        # - basis_strategy: Options include 'residuals', 'jacobian', 'reactions'
        # - solving_technique: Options include 'normal_equations', 'qr_decomposition'
        rom_bns_settings = """{
            "train_petrov_galerkin": false,
            "basis_strategy": "residuals",
            "include_phi": false,
            "svd_truncation_tolerance": 1e-8,
            "solving_technique": "normal_equations",
            "monotonicity_preserving": false
        }"""
        return rom_bns_settings



    def kmeans_method(self, test_data, n_clusters, mu, cluster_list):
        if n_clusters > len(mu): 
            n_clusters = len(mu)
            print("Number of clusters corrected to ", len(mu))
        kmeans_object = KMeans(n_clusters=n_clusters, random_state=0).fit(test_data)
        cluster_list.append(kmeans_object.labels_)
        fig = plt.figure()
        fig.set_figwidth(12.0)
        fig.set_figheight(8.0)
        self.centroids_to_plot = np.zeros([n_clusters,2])
        for j in range(n_clusters):
            fig = plt.scatter(mu[:, 1][kmeans_object.labels_==j], mu[:, 0][kmeans_object.labels_==j], label="Cluster "+str(j))
            self.centroids_to_plot[j,:] = np.mean(mu[:, 1][kmeans_object.labels_==j]), np.mean(mu[:, 0][kmeans_object.labels_==j])
            fig = plt.scatter(self.centroids_to_plot[j,0], self.centroids_to_plot[j,1], c='k',marker="s", s= 150)
            fig = plt.text(self.centroids_to_plot[j,0], self.centroids_to_plot[j,1]-0.05, "Centroid C"+str(j), ha='center')
            #save the centroid data (snapshotmatrix)
            # np.save("k-means_initial_solution_cluster_"+str(j)+"_"+str(self.centroids_to_plot[j,1])+"_"+str(self.centroids_to_plot[j,0]),
            #         kmeans_object.cluster_centers_.T[:,j])
        # fig = plt.scatter(self.centroids_to_plot[:,0], self.centroids_to_plot[:,1], c='k',marker="s", s= 150)
        fig = plt.title("k-means clustering")
        fig = plt.ylabel('Alpha')
        fig = plt.xlabel('Mach')
        # fig = plt.minorticks_on()
        # fig = plt.grid(which='both')
        fig = plt.legend(bbox_to_anchor=(1.004, 0.9, 1.0, 0.102), loc='upper left', borderaxespad=0.)
        # fig = plt.show()
        fig = plt.savefig('k-means clustering.pdf',bbox_inches='tight' )
        fig = plt.close('all')
        # To obtein the number of the cluster
        # print("Cluster for '1.448, 0.732.npy': " + str(kmeans_object.predict(np.load("1.448, 0.732.npy").T)))
        return n_clusters


    def E_p(self, u, c):
        """
        c: direction vector onto which to project
        u: vector or colection of column vectors to project onto the direction of c
        """
        c = c.reshape(-1,1)
        if len(u.shape)==1:
            u = u.reshape(-1,1)
        projection_of_u_onto_c = ((c@c.T) / (c.T@c)) @ u
        projection_error = np.linalg.norm(u - projection_of_u_onto_c, axis=0) / np.linalg.norm(u,axis=0)
        return projection_error

    def PEBL(self, Snapshots, bisection_tolerance,  POD_tolerance):
        #stage 1, generation of bisection tree with accuracy 'bisection_tolerance'
        max_index = np.argmax( np.linalg.norm(Snapshots, axis=0) )
        first_snapshot = Snapshots[:,max_index]
        Tree = Node([first_snapshot, Snapshots])
        bisect_flag = True
        while bisect_flag == True:
            bisect_flag = False
            for leaf in Tree.leaves():
                errors = self.E_p(leaf.val[1], leaf.val[0])
                max_error = max(errors)
                print(max_error)
                if max_error > bisection_tolerance:
                    bisect_flag = True
                    #find next anchor point
                    max_index = np.argmax(errors)
                    c_new = leaf.val[1][:,max_index]
                    new_errors = self.E_p(leaf.val[1], c_new)
                    indexes_left = np.where( errors <= new_errors)
                    indexes_right = np.where( errors > new_errors)
                    #divide the snapshots among the two children
                    leaf.left =  Node([leaf.val[0], leaf.val[1][:,indexes_left[0]]])
                    leaf.right = Node([c_new, leaf.val[1][:,indexes_right[0]]])
                    leaf.val[1] = None
        #stage 2, generation of local POD bases'
        for leaf in Tree.leaves():
            Phi_i = self.ObtainBasis(leaf.val[1], POD_tolerance)
            leaf.val.append(Phi_i)
        return Tree


    def ObtainBasis(self, Snapshots, truncation_tolerance=0):
            U,_,_= self.truncated_svd(Snapshots,truncation_tolerance)
            return U

    def truncated_svd(self, Matrix, epsilon=0):
        M,N=np.shape(Matrix)
        dimMATRIX = max(M,N)
        U, s, V = np.linalg.svd(Matrix, full_matrices=False) #U --> M xN, V --> N x N
        V = V.T
        tol = dimMATRIX*np.finfo(float).eps*max(s)/2
        R = np.sum(s > tol)  # Definition of numerical rank
        if epsilon == 0:
            K = R
        else:
            SingVsq = np.multiply(s,s)
            SingVsq.sort()
            normEf2 = np.sqrt(np.cumsum(SingVsq))
            epsilon = epsilon*normEf2[-1] #relative tolerance
            T = (sum(normEf2<epsilon))
            K = len(s)-T
        K = min(R,K)
        return U[:, :K], s[:K], V[:, :K]

    def pebl_method(self, test_data,  bisection_tolerance,  POD_tolerance, cluster_list=False, list=False):
        Tree = self.PEBL(test_data.T,  bisection_tolerance=bisection_tolerance,  POD_tolerance=POD_tolerance)
        fig = plt.figure()
        fig.set_figwidth(12.0)
        fig.set_figheight(8.0)
        print("leaves: ", len(Tree.leaves()))
        self.centroids_to_plot = np.zeros([len(Tree.leaves()),2])
        clusters = np.zeros(len(list))
        leaf_counter = 0
        for leaf in Tree.leaves():
            #save the centroid data (snapshotmatrix)
            # np.save("initial_solution_cluster_"+str(leaf_counter),leaf.val[0])
            print("Number of cases in leaf "+str(leaf_counter)+": ", len(leaf.val[1][1]))
            angle = np.zeros(len(leaf.val[1][1]))
            mach  = np.zeros(len(leaf.val[1][1]))
            for j in range(len(leaf.val[1][1])): 
                for i in range(test_data.T.shape[1]):
                    if np.array_equal(test_data.T[:,i], leaf.val[1][:,j]):
                        angle[j]    = list[i][0]
                        mach[j]     = list[i][1]
                        clusters[i] = np.int0(leaf_counter)
                        break
            fig = plt.scatter(mach, angle, label="Cluster "+str(leaf_counter))
            self.centroids_to_plot[leaf_counter,:] = np.mean(mach), np.mean(angle)
            fig = plt.scatter(self.centroids_to_plot[leaf_counter,0], self.centroids_to_plot[leaf_counter,1], c='k',marker="s", s= 150)
            fig = plt.text(self.centroids_to_plot[leaf_counter,0], self.centroids_to_plot[leaf_counter,1]-0.05, "Centroid C"+str(leaf_counter), ha='center')
            #save the centroid data (snapshotmatrix)
            # np.save("pebl_initial_solution_cluster_"+str(leaf_counter)+"_"+str(np.mean(angle))+"_"+str(np.mean(mach)),leaf.val[0])
            leaf_counter+=1
        # fig = plt.scatter(self.centroids_to_plot[:,0], self.centroids_to_plot[:,1], c='k',marker="s", s= 150)
        cluster_list.append(clusters)
        fig = plt.title("PEBL clustering")
        fig = plt.ylabel('Alpha')
        fig = plt.xlabel('Mach')
        # fig = plt.minorticks_on()
        # fig = plt.grid(which='both')
        fig = plt.legend(bbox_to_anchor=(1.004, 0.9, 1.0, 0.102), loc='upper left', borderaxespad=0.)
        fig = plt.savefig('PEBL clustering.pdf',bbox_inches='tight' )
        # fig = plt.show()
        fig = plt.close('all')
        return len(Tree.leaves())

    def PlotClusterListWithRomModes(self,n_clusters,mu_list,cluster_list,rom_modes_list,clustering_method):
        fig = plt.figure()
        fig.set_figwidth(15.0)
        fig.set_figheight(10.0)
        for i in range(n_clusters):
            mu_aux = [mu for mu, cluster in zip(mu_list, cluster_list[:][0]) if cluster == i]
            mu_train_a = np.zeros(len(mu_aux))
            mu_train_m = np.zeros(len(mu_aux))
            for j in range(len(mu_aux)):
                mu_train_a[j] = mu_aux[j][0]
                mu_train_m[j] = mu_aux[j][1]
            fig = plt.scatter(mu_train_m, mu_train_a, label="Cluster " + str(i) + ": Rom modes: " + str(np.int0(rom_modes_list[i])))
            fig = plt.scatter(self.centroids_to_plot[i,0], self.centroids_to_plot[i,1], c='k',marker="s", s= 150)
            fig = plt.text(self.centroids_to_plot[i,0], self.centroids_to_plot[i,1]-0.05, "Centroid C"+str(i), ha='center')
        # fig = plt.scatter(self.centroids_to_plot[:,0], self.centroids_to_plot[:,1], c='k',marker="s", s= 150)
        fig = plt.title(clustering_method+" clustering")
        fig = plt.ylabel('Alpha')
        fig = plt.xlabel('Mach')
        # fig = plt.minorticks_on()
        # fig = plt.grid(which='both')
        fig = plt.legend(bbox_to_anchor=(1.004, 0.9, 1.0, 0.102), loc='upper left', borderaxespad=0.)
        # fig = plt.show()
        fig = plt.savefig(clustering_method+' clustering.pdf',bbox_inches='tight' )
        fig = plt.close('all')

    def pebl_method_full_snapshots_matrix(self, test_data, simple_matrix,  bisection_tolerance,  POD_tolerance, cluster_list=False, list=False):
        Tree = self.PEBL(test_data.T,  bisection_tolerance=bisection_tolerance,  POD_tolerance=POD_tolerance)
        print("leaves: ", len(Tree.leaves()))
        self.centroids_to_plot = np.zeros([len(Tree.leaves()),2])
        clusters = np.zeros(len(list))
        leaf_counter = 0
        for leaf in Tree.leaves():
            print("Number of cases in leaf "+str(leaf_counter)+": ", len(leaf.val[1][1]))
            for j in range(len(leaf.val[1][1])): 
                for i in range(simple_matrix.T.shape[1]):
                    if np.array_equal(simple_matrix.T[:,i], leaf.val[1][:,j]):
                        clusters[i] = leaf_counter
                        break
            leaf_counter+=1
        cluster_list.append(clusters)
        return len(Tree.leaves())
    
    def GetClusterId(self, centroids_list, mu):
        return np.argmin(np.sqrt(np.sum((centroids_list - mu) ** 2, axis=1)))

class Node:
    def __init__(self, val):
        self.left = None
        self.right = None
        self.val = val

    def leaves(self):
        current_nodes = [self]
        leaves = []

        while len(current_nodes) > 0:
            next_nodes = []
            for node in current_nodes:
                if node.left is None and node.right is None:
                    leaves.append(node)
                    continue
                if node.left is not None:
                    next_nodes.append(node.left)
                if node.right is not None:
                    next_nodes.append(node.right)
            current_nodes = next_nodes
        return leaves