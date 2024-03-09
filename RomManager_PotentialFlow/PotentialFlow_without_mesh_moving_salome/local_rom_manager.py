import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis import PotentialFlowAnalysis
from KratosMultiphysics.RomApplication.rom_testing_utilities import SetUpSimulationInstance
from KratosMultiphysics.RomApplication.calculate_rom_basis_output_process import CalculateRomBasisOutputProcess
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition
import numpy as np
import os
import importlib
import json
from pathlib import Path



class LocalRomManager(object):

    def __init__(self,project_parameters_name, general_rom_manager_parameters, CustomizeSimulation, UpdateProjectParameters):
        self.project_parameters_name = project_parameters_name
        self.general_rom_manager_parameters = general_rom_manager_parameters
        self.CustomizeSimulation = CustomizeSimulation
        self.UpdateProjectParameters = UpdateProjectParameters
        self.ROMvsFOM_train     = 0.0
        self.ROMvsHROM_train    = 0.0
        self.ROMvsFOM_test      = 0.0
        self.ROMvsHROM_test     = 0.0


    def Fit(self, mu_train=[None], use_uncorrected_solutions=False):
        mu_train_errors = []
        if (len(mu_train) > 0):
            chosen_projection_strategy = self.general_rom_manager_parameters["projection_strategy"].GetString()
            training_stages = self.general_rom_manager_parameters["rom_stages_to_train"].GetStringArray()
            #######################
            ######  Galerkin ######
            if chosen_projection_strategy == "galerkin":
                if any(item == "ROM" for item in training_stages):
                    fom_snapshots = self.__LaunchTrainROM(mu_train, use_uncorrected_solutions)
                    self._ChangeRomFlags(simulation_to_run = "GalerkinROM")
                    rom_snapshots = self.__LaunchROM(mu_train)
                    self.ROMvsFOM_train = np.linalg.norm(fom_snapshots - rom_snapshots)*100/ np.linalg.norm(fom_snapshots)
                    for i in range(len(mu_train)):
                        error = np.linalg.norm(fom_snapshots[:,i] - rom_snapshots[:,i])*100/ np.linalg.norm(fom_snapshots[:,i])
                        mu_train_errors.append([mu_train[i][0],mu_train[i][1],error])

                if any(item == "HROM" for item in training_stages):
                    #FIXME there will be an error if we only train HROM, but not ROM
                    self._ChangeRomFlags(simulation_to_run = "trainHROMGalerkin")
                    self.__LaunchTrainHROM(mu_train)
                    self._ChangeRomFlags(simulation_to_run = "runHROMGalerkin")
                    hrom_snapshots = self.__LaunchHROM(mu_train)
                    self.ROMvsHROM_train = np.linalg.norm(rom_snapshots - hrom_snapshots)*100/ np.linalg.norm(rom_snapshots)
                    # for i in range(len(mu_train)):
                    #     error = np.linalg.norm(rom_snapshots[:,i] - hrom_snapshots[:,i])*100/ np.linalg.norm(rom_snapshots[:,i])
                    #     mu_train_errors.append([mu_train[i][0],mu_train[i][1],error])
            #######################

            #######################################
            ##  Least-Squares Petrov Galerkin   ###
            elif chosen_projection_strategy == "lspg":
                if any(item == "ROM" for item in training_stages):
                    fom_snapshots = self.__LaunchTrainROM(mu_train, use_uncorrected_solutions)
                    self._ChangeRomFlags(simulation_to_run = "lspg")
                    rom_snapshots = self.__LaunchROM(mu_train)
                    self.ROMvsFOM_train = np.linalg.norm(fom_snapshots - rom_snapshots)*100/ np.linalg.norm(fom_snapshots)
                    for i in range(len(mu_train)):
                        error = np.linalg.norm(fom_snapshots[:,i] - rom_snapshots[:,i])*100/ np.linalg.norm(fom_snapshots[:,i])
                        mu_train_errors.append([mu_train[i][0],mu_train[i][1],error])
                if any(item == "HROM" for item in training_stages):
                    # Change the flags to train the HROM for LSPG
                    self._ChangeRomFlags(simulation_to_run = "trainHROMLSPG")
                    self.__LaunchTrainHROM(mu_train)

                    # Change the flags to run the HROM for LSPG
                    self._ChangeRomFlags(simulation_to_run = "runHROMLSPG")
                    hrom_snapshots = self.__LaunchHROM(mu_train)
                    self.ROMvsHROM_train = np.linalg.norm(rom_snapshots - hrom_snapshots)*100/ np.linalg.norm(rom_snapshots)
                    # for i in range(len(mu_train)):
                    #     error = np.linalg.norm(rom_snapshots[:,i] - hrom_snapshots[:,i])*100/ np.linalg.norm(rom_snapshots[:,i])
                    #     mu_train_errors.append([mu_train[i][0],mu_train[i][1],error])
                    #######################################

            ##########################
            ###  Petrov Galerkin   ###
            elif chosen_projection_strategy == "petrov_galerkin":
                if any(item == "ROM" for item in training_stages):
                    fom_snapshots = self.__LaunchTrainROM(mu_train, use_uncorrected_solutions)
                    self._ChangeRomFlags(simulation_to_run = "TrainPG")
                    self.__LaunchTrainPG(mu_train, use_uncorrected_solutions)
                    self._ChangeRomFlags(simulation_to_run = "PG")
                    rom_snapshots = self.__LaunchROM(mu_train)
                    self.ROMvsFOM_train = np.linalg.norm(fom_snapshots - rom_snapshots)*100/ np.linalg.norm(fom_snapshots)
                    for i in range(len(mu_train)):
                        error = np.linalg.norm(fom_snapshots[:,i] - rom_snapshots[:,i])*100/ np.linalg.norm(fom_snapshots[:,i])
                        mu_train_errors.append([mu_train[i][0],mu_train[i][1],error])
                if any(item == "HROM" for item in training_stages):
                    #FIXME there will be an error if we only train HROM, but not ROM
                    self._ChangeRomFlags(simulation_to_run = "trainHROMPetrovGalerkin")
                    self.__LaunchTrainHROM(mu_train)
                    self._ChangeRomFlags(simulation_to_run = "runHROMPetrovGalerkin")
                    hrom_snapshots = self.__LaunchHROM(mu_train)
                    self.ROMvsHROM_train = np.linalg.norm(rom_snapshots - hrom_snapshots)*100/ np.linalg.norm(rom_snapshots)
                    # for i in range(len(mu_train)):
                    #     error = np.linalg.norm(rom_snapshots[:,i] - hrom_snapshots[:,i])*100/ np.linalg.norm(rom_snapshots[:,i])
                    #     mu_train_errors.append([mu_train[i][0],mu_train[i][1],error])
            ##########################
            else:
                err_msg = f'Provided projection strategy {chosen_projection_strategy} is not supported. Available options are \'galerkin\', \'lspg\' and \'petrov_galerkin\'.'
                raise Exception(err_msg)
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
                    self._ChangeRomFlags(simulation_to_run = "GalerkinROM")
                    rom_snapshots = self.__LaunchTestROM(mu_test)
                    self.ROMvsFOM_test = np.linalg.norm(fom_snapshots - rom_snapshots)*100/ np.linalg.norm(fom_snapshots)
                    for i in range(len(mu_test)):
                        error = np.linalg.norm(fom_snapshots[:,i] - rom_snapshots[:,i])*100/ np.linalg.norm(fom_snapshots[:,i])
                        mu_test_errors.append([mu_test[i][0],mu_test[i][1],error])

                if any(item == "HROM" for item in testing_stages):
                    #FIXME there will be an error if we only test HROM, but not ROM
                    self._ChangeRomFlags(simulation_to_run = "runHROMGalerkin")
                    hrom_snapshots = self.__LaunchTestHROM(mu_test)
                    self.ROMvsHROM_test = np.linalg.norm(rom_snapshots - hrom_snapshots)*100/ np.linalg.norm(rom_snapshots)
                    # for i in range(len(mu_test)):
                    #     error = np.linalg.norm(rom_snapshots[:,i] - hrom_snapshots[:,i])*100/ np.linalg.norm(rom_snapshots[:,i])
                    #     mu_test_errors.append([mu_test[i][0],mu_test[i][1],error])


            #######################################
            ##  Least-Squares Petrov Galerkin   ###
            elif chosen_projection_strategy == "lspg":
                if any(item == "ROM" for item in testing_stages):
                    fom_snapshots = self.__LaunchTestFOM(mu_test)
                    self._ChangeRomFlags(simulation_to_run = "lspg")
                    rom_snapshots = self.__LaunchTestROM(mu_test)
                    self.ROMvsFOM_test = np.linalg.norm(fom_snapshots - rom_snapshots)*100/ np.linalg.norm(fom_snapshots)
                    for i in range(len(mu_test)):
                        error = np.linalg.norm(fom_snapshots[:,i] - rom_snapshots[:,i])*100/ np.linalg.norm(fom_snapshots[:,i])
                        mu_test_errors.append([mu_test[i][0],mu_test[i][1],error])
                if any(item == "HROM" for item in testing_stages):
                    self._ChangeRomFlags(simulation_to_run = "runHROMLSPG")
                    hrom_snapshots = self.__LaunchTestHROM(mu_test)
                    self.ROMvsHROM_test = np.linalg.norm(rom_snapshots - hrom_snapshots)*100/ np.linalg.norm(rom_snapshots)
                    # for i in range(len(mu_test)):
                    #     error = np.linalg.norm(rom_snapshots[:,i] - hrom_snapshots[:,i])*100/ np.linalg.norm(rom_snapshots[:,i])
                    #     mu_test_errors.append([mu_test[i][0],mu_test[i][1],error])
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
                        mu_test_errors.append([mu_test[i][0],mu_test[i][1],error])
                if any(item == "HROM" for item in testing_stages):
                    #FIXME there will be an error if we only train HROM, but not ROM
                    self._ChangeRomFlags(simulation_to_run = "runHROMPetrovGalerkin")
                    hrom_snapshots = self.__LaunchTestHROM(mu_test)
                    self.ROMvsHROM_test = np.linalg.norm(rom_snapshots - hrom_snapshots)*100/ np.linalg.norm(rom_snapshots)
                    # for i in range(len(mu_test)):
                    #     error = np.linalg.norm(rom_snapshots[:,i] - hrom_snapshots[:,i])*100/ np.linalg.norm(rom_snapshots[:,i])
                    #     mu_test_errors.append([mu_test[i][0],mu_test[i][1],error])
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
        if any(item == "ROM" for item in testing_stages) and self.ROMvsFOM_test != 0:
            print("approximation error in test set FOM vs ROM: ", self.ROMvsFOM_test," %")
        if any(item == "HROM" for item in testing_stages) and self.ROMvsHROM_test != 0:
            print("approximation error in test set ROM vs HROM: ",  self.ROMvsHROM_test," %")


    def __LaunchTrainROM(self, mu_train, use_uncorrected_solutions):
        """
        This method should be parallel capable
        """     
        if not os.path.exists("DataBase"):
            os.mkdir("DataBase")
            os.mkdir("DataBase/Results")
            os.mkdir("DataBase/Data")
            os.mkdir("DataBase/UncorrectedSolutions")
            os.mkdir("DataBase/UncorrectedSolutions/Data")
            os.mkdir("DataBase/UncorrectedSolutions/Results")
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        SnapshotsMatrix = []
        error_SnapshotsMatrix = []
        for Id, mu in enumerate(mu_train):
            aux_snapshotmatrix = []
            parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
            parameters_copy = self._AddBasisCreationToProjectParameters(parameters_copy) #TODO stop using the RomBasisOutputProcess to store the snapshots. Use instead the upcoming build-in function
            parameters_copy = self._StoreResultsByName(parameters_copy,'FOM_Fit',mu,Id)
            if os.path.exists("DataBase/Results/" + str(mu[0]) + ", " + str(mu[1])):
                parameters_copy["output_processes"].RemoveValue("vtk_output")
            else:
                parameters_copy["output_processes"]["vtk_output"][0]["Parameters"]["output_path"].SetString("DataBase/Results/"+ str(mu[0]) + ", " + str(mu[1]))
            model = KratosMultiphysics.Model()
            analysis_stage_class = self._GetAnalysisStageClass(parameters_copy)
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy)
            case_name = "DataBase/"+ str(mu[0]) + ", " + str(mu[1]) +".npy"
            if os.path.exists(case_name):
                aux_snapshotmatrix = np.load(case_name)
                simulation.Initialize()
                for process in simulation._GetListOfOutputProcesses():
                    if isinstance(process, CalculateRomBasisOutputProcess):
                        BasisOutputProcess = process
            else:
                simulation.Run()
                for process in simulation._GetListOfOutputProcesses():
                    if isinstance(process, CalculateRomBasisOutputProcess):
                        BasisOutputProcess = process
                aux_snapshotmatrix = BasisOutputProcess._GetSnapshotsMatrix()
                np.save(case_name, aux_snapshotmatrix)
            if use_uncorrected_solutions:
                Uncorrected_case_name = "DataBase/UncorrectedSolutions/"+ str(mu[0]) + ", " + str(mu[1]) +".npy"
                if os.path.exists(Uncorrected_case_name):
                    aux_Uncorrected_snapshotmatrix = np.load(Uncorrected_case_name)
                    SnapshotsMatrix.append(aux_Uncorrected_snapshotmatrix)
                else:
                    aux_parameters = parameters_copy.Clone()
                    aux_parameters["solver_settings"]["scheme_settings"]["update_critical_mach"].SetDouble(0.85)
                    aux_parameters["solver_settings"]["scheme_settings"]["update_upwind_factor_constant"].SetDouble(2.0)
                    aux_parameters["solver_settings"]["scheme_settings"]["update_transonic_tolerance"].SetDouble(1e-30)
                    if os.path.exists("DataBase/UncorrectedSolutions/Results/" + str(mu[0]) + ", " + str(mu[1])):
                        aux_parameters["output_processes"].RemoveValue("vtk_output")
                    else:
                        aux_parameters["output_processes"]["vtk_output"][0]["Parameters"]["output_path"].SetString("DataBase/UncorrectedSolutions/Results/"+ str(mu[0]) + ", " + str(mu[1]))
                    aux_model = KratosMultiphysics.Model()
                    aux_simulation = PotentialFlowAnalysis(aux_model,aux_parameters)
                    aux_simulation.Run()
                    for process in aux_simulation._GetListOfOutputProcesses():
                        if isinstance(process, CalculateRomBasisOutputProcess):
                            aux_BasisOutputProcess = process
                    aux_Uncorrected_snapshotmatrix = aux_BasisOutputProcess._GetSnapshotsMatrix()
                    np.save(Uncorrected_case_name, aux_Uncorrected_snapshotmatrix)
                    SnapshotsMatrix.append(aux_Uncorrected_snapshotmatrix)

                    skin_data_name = "DataBase/UncorrectedSolutions/Data/"+ str(mu[0]) + ", " + str(mu[1]) +".dat"
                    fout = open(skin_data_name,'w')
                    modelpart = aux_model["MainModelPart.Body2D_Body"]
                    for node in modelpart.Nodes:
                        x = node.X ; y = node.Y ; z = node.Z
                        cp = node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
                        fout.write("%s %s %s %s\n" %(x,y,z,cp))
                    fout.close()
                    fout=open("DataBase/UncorrectedSolutions/Data/full_" + str(mu[0]) + ", " + str(mu[1]) + ".dat",'w')
                    modelpart = aux_model["MainModelPart"]
                    for node in modelpart.Nodes:
                        id_node = node.Id
                        velocity_potential = node.GetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL)
                        auxiliary_velocity_potential = node.GetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL)
                        fout.write("%s %s %s\n" %(id_node, velocity_potential, auxiliary_velocity_potential))
                    fout.close()
            SnapshotsMatrix.append(aux_snapshotmatrix) #TODO add a CustomMethod() as a standard method in the Analysis Stage to retrive some solution
            error_SnapshotsMatrix.append(aux_snapshotmatrix)
        SnapshotsMatrix = np.block(SnapshotsMatrix)
        error_SnapshotsMatrix = np.block(error_SnapshotsMatrix)
        BasisOutputProcess._PrintRomBasis(SnapshotsMatrix) #Calling the RomOutput Process for creating the RomParameter.json

        return error_SnapshotsMatrix


    def __LaunchROM(self, mu_train):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        SnapshotsMatrix = []
        for Id, mu in enumerate(mu_train):
            parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
            parameters_copy = self._AddBasisCreationToProjectParameters(parameters_copy)  #TODO stop using the RomBasisOutputProcess to store the snapshots. Use instead the upcoming build-in function
            parameters_copy = self._StoreResultsByName(parameters_copy,'ROM_Fit',mu,Id)
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


    def __LaunchTrainPG(self, mu_train, use_uncorrected_solutions):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        PetrovGalerkinTrainMatrix = []
        for Id, mu in enumerate(mu_train):
            parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
            parameters_copy = self._AddBasisCreationToProjectParameters(parameters_copy) 
            # parameters_copy = self._StoreResultsByName(parameters_copy,'FOM_Fit',mu,Id)
            parameters_copy = self._StoreNoResults(parameters_copy)
            model = KratosMultiphysics.Model()
            analysis_stage_class = type(SetUpSimulationInstance(model, parameters_copy))
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy)
            simulation.Run()
            PetrovGalerkinTrainMatrix.append(simulation.GetPetrovGalerkinTrainUtility()._GetSnapshotsMatrix()) #TODO is the best way of extracting the Projected Residuals calling the HROM residuals utility?
            
            if use_uncorrected_solutions:
                parameters_copy["solver_settings"]["scheme_settings"]["update_critical_mach"].SetDouble(0.85)
                parameters_copy["solver_settings"]["scheme_settings"]["update_upwind_factor_constant"].SetDouble(2.0)
                parameters_copy["solver_settings"]["scheme_settings"]["update_transonic_tolerance"].SetDouble(1e-30)
                aux_model = KratosMultiphysics.Model()
                analysis_stage_class = type(SetUpSimulationInstance(aux_model, parameters_copy))
                aux_simulation = self.CustomizeSimulation(analysis_stage_class,aux_model,parameters_copy)
                aux_simulation.Run()
                PetrovGalerkinTrainMatrix.append(simulation.GetPetrovGalerkinTrainUtility()._GetSnapshotsMatrix())

        simulation.GetPetrovGalerkinTrainUtility().CalculateAndSaveBasis(np.block(PetrovGalerkinTrainMatrix))


    def __LaunchTrainHROM(self, mu_train):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        RedidualsSnapshotsMatrix = []
        for mu in mu_train:
            parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
            parameters_copy = self._AddBasisCreationToProjectParameters(parameters_copy)
            parameters_copy = self._StoreNoResults(parameters_copy)
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


    def __LaunchHROM(self, mu_train):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        SnapshotsMatrix = []
        for Id, mu in enumerate(mu_train):
            parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
            parameters_copy = self._AddBasisCreationToProjectParameters(parameters_copy)
            parameters_copy = self._StoreResultsByName(parameters_copy,'HROM_Fit',mu,Id)
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
        if not os.path.exists("DataBase"):
            os.mkdir("DataBase")
            os.mkdir("DataBase/Results")
            os.mkdir("DataBase/Data")
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        SnapshotsMatrix = []
        for Id, mu in enumerate(mu_test):
            aux_snapshotmatrix = []
            parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
            parameters_copy = self._AddBasisCreationToProjectParameters(parameters_copy) #TODO stop using the RomBasisOutputProcess to store the snapshots. Use instead the upcoming build-in function
            parameters_copy = self._StoreResultsByName(parameters_copy,'FOM_Test',mu,Id)
            if os.path.exists("DataBase/Results/" + str(mu[0]) + ", " + str(mu[1])):
                parameters_copy["output_processes"].RemoveValue("vtk_output")
            else:
                parameters_copy["output_processes"]["vtk_output"][0]["Parameters"]["output_path"].SetString("DataBase/Results/"+ str(mu[0]) + ", " + str(mu[1]))
            model = KratosMultiphysics.Model()
            analysis_stage_class = self._GetAnalysisStageClass(parameters_copy)
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy)
            case_name = "DataBase/"+ str(mu[0]) + ", " + str(mu[1]) +".npy"
            if os.path.exists(case_name):
                aux_snapshotmatrix = np.load(case_name)
                simulation.Initialize()
                for process in simulation._GetListOfOutputProcesses():
                    if isinstance(process, CalculateRomBasisOutputProcess):
                        BasisOutputProcess = process
            else:
                simulation.Run()
                for process in simulation._GetListOfOutputProcesses():
                    if isinstance(process, CalculateRomBasisOutputProcess):
                        BasisOutputProcess = process
                aux_snapshotmatrix = BasisOutputProcess._GetSnapshotsMatrix()
                np.save(case_name, aux_snapshotmatrix)
            SnapshotsMatrix.append(aux_snapshotmatrix) #TODO add a CustomMethod() as a standard method in the Analysis Stage to retrive some solution
        SnapshotsMatrix = np.block(SnapshotsMatrix)

        return SnapshotsMatrix


    def __LaunchTestROM(self, mu_test):
        """
        This method should be parallel capable
        """
        SnapshotsMatrix = []
        for Id, mu in enumerate(mu_test):
            with open(self.project_parameters_name,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
            parameters_copy = self._AddBasisCreationToProjectParameters(parameters_copy)  #TODO stop using the RomBasisOutputProcess to store the snapshots. Use instead the upcoming build-in function
            parameters_copy = self._StoreResultsByName(parameters_copy,'ROM_Test',mu,Id)
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