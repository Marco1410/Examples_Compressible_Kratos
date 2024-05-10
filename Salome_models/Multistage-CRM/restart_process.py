# Importing the Kratos Library
import KratosMultiphysics

# Other imports
import os
from pathlib import Path

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return RestartProcess(model, settings["Parameters"])

# All the processes python should be derived from "Process"
class RestartProcess(KratosMultiphysics.Process):
    """This process allow save and load restart files each step or time of solving selected."
    
    Public member variables:
    model    -- The container of the different model parts.
    settings -- Kratos parameters containing process settings.
    """

    def __init__(self, model, settings):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing process settings.
        """
        
        # calling the baseclass constructor
        KratosMultiphysics.Process.__init__(self) 

        default_settings = KratosMultiphysics.Parameters("""{
                                "save_restart"                   : false,
                                "load_restart"                   : false,
                                "clean_previus_files"            : false,
                                "model_part"                     : "FluidModelPart.fluid_computational_model_part",
                                "input_filename"                 : "Restart",
                                "input_output_path"              : "restart_folder",
                                "echo_level"                     : 1,
                                "serializer_trace"               : "no_trace",
                                "restart_load_file_label"        : "",
                                "load_restart_files_from_folder" : true,
                                "restart_save_frequency"         : 0.0,
                                "restart_control_type"           : "time",
                                "save_restart_files_in_folder"   : true,
                                "max_files_to_keep"              : -1
        }""")

        settings.ValidateAndAssignDefaults(default_settings)
        self.save_restart = settings["save_restart"].GetBool()
        self.load_restart = settings["load_restart"].GetBool()
        self.model_part_name = settings["model_part"].GetString()
        self.clean_previus_files = settings["clean_previus_files"].GetBool()
        settings.RemoveValue("model_part")
        settings.RemoveValue("save_restart")
        settings.RemoveValue("load_restart")
        settings.RemoveValue("clean_previus_files")
        self.max_files_to_keep = settings["max_files_to_keep"].GetInt()
        self.restart_control_type = settings["restart_control_type"].GetString()
        self.save_restart_files_in_folder = settings["save_restart_files_in_folder"].GetBool()    
        self.raw_path, self.raw_file_name = os.path.split(settings["input_filename"].GetString())
        self.raw_path = os.path.join(os.getcwd(), self.raw_path)
        if settings["input_output_path"].GetString() == '':
            self.input_output_path = self.raw_file_name + "__restart_files"
        else:
            self.input_output_path = settings["input_output_path"].GetString()
        
        # saving the modelpart
        self.model_part = model[self.model_part_name]
        
        if self.model_part.IsDistributed(): # mpi-execution
            from KratosMultiphysics.mpi.distributed_restart_utility import DistributedRestartUtility as RestartUtility
        else:
            from KratosMultiphysics.restart_utility import RestartUtility

        self.restart_utility = RestartUtility(self.model_part, settings)
     
    def ExecuteBeforeSolutionLoop(self):
        if self.save_restart:
            if self.clean_previus_files:
                self._CleanPreviusFiles()
            self.restart_files_number = 0
            self.restart_files = {}  
        if self.load_restart:
            self.restart_utility.LoadRestart()
        
    def ExecuteFinalizeSolutionStep(self):
        if self.save_restart:
            if self.restart_utility.IsRestartOutputStep():
                self.restart_utility.SaveRestart()
                self.restart_files_number += 1
                if self.restart_control_type == "time":
                    time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
                    self.control_label = self.__GetPrettyTime(time)
                elif self.restart_control_type == "step":
                    self.control_label = self.model_part.ProcessInfo[KratosMultiphysics.STEP]
                file_name = self.__GetFileNameSave(self.control_label)
                self.restart_files[self.control_label] = file_name + ".rest"
                self._ClearObsoleteRestartFiles()

    def _ClearObsoleteRestartFiles(self):
        """Delete restart files that are no longer needed."""
        if self.max_files_to_keep > -1:

            number_of_obsolete_files = self.restart_files_number - self.max_files_to_keep
            restart_file_keys = sorted(self.restart_files)
            
            for i in range(number_of_obsolete_files):
                i_key = restart_file_keys[i]
                file_path = os.path.join(self.__GetFolderPathSave(), self.restart_files[i_key])
                KratosMultiphysics.kratos_utilities.DeleteFileIfExisting(file_path)          

    def _CleanPreviusFiles(self):
        self._CreateOutputFolder()
        path = self.__GetFolderPathSave() + "/"
        filelist = [f for f in os.listdir(path) if f.endswith(".rest")]
        for f in filelist:
            try:
                KratosMultiphysics.kratos_utilities.DeleteFileIfExisting(os.path.join(self.__GetFolderPathSave(), f))
            except OSError:
                pass
                
    def _CreateOutputFolder(self):
        if self.save_restart_files_in_folder:
            KratosMultiphysics.FilesystemExtensions.MPISafeCreateDirectories(self.__GetFolderPathSave())
            
    def __GetFolderPathSave(self):
        if self.save_restart_files_in_folder:
            return os.path.join(self.raw_path, self.input_output_path)
        else:
            return self.raw_path
                
    def _GetFileLabelSave(self, file_label):
        return str(file_label) 
        
    def __GetFileNameSave(self, file_label):
        restart_file_name = self.raw_file_name + '_' + self._GetFileLabelSave(file_label)
        return os.path.join(self.__GetFolderPathSave(), restart_file_name)
          
    def __GetPrettyTime(self, time):
        """This functions reduces the digits of a number to a relevant precision."""
        pretty_time = "{0:.12g}".format(time)
        pretty_time = float(pretty_time)
        return pretty_time
