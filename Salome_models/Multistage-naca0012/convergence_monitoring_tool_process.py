# Importing the Kratos Library
import KratosMultiphysics

# Other imports
import sys
import math

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ConvergenceMonitoringToolProcess(model, settings["Parameters"])

# All the processes python should be derived from "Process"
class ConvergenceMonitoringToolProcess(KratosMultiphysics.Process):
    """This process allow convergence monitoring in each step of solving."
    
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
                                "print_error_estimation" : false,
                                "print_warnings"         : false,
                                "tolerance"              : 1e-12,
                                "model_part"             : "FluidModelPart.fluid_computational_model_part"
        }""")

        settings.ValidateAndAssignDefaults(default_settings)
        self.print_warnings = settings["print_warnings"].GetBool()
        self.print_error_estimation = settings["print_error_estimation"].GetBool()
        self.model_part_name = settings["model_part"].GetString()
        self.tolerance = settings["tolerance"].GetDouble()
                
        # saving the modelpart
        self.model_part = model[self.model_part_name]
    
    def ExecuteBeforeSolutionLoop(self):
        KratosMultiphysics.Logger.Print("::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
            
    def ExecuteFinalizeSolutionStep(self):
        
        step = self.model_part.ProcessInfo.GetValue(KratosMultiphysics.STEP)
        
        if  step > 1 :

            e_density     = 0  
            e_momentum    = 0
            e_totalEnergy = 0
            n_density     = 0  
            n_momentum    = 0
            n_totalEnergy = 0

            for node in self.model_part.Nodes:
                e_density     = e_density + (node.GetSolutionStepValue(KratosMultiphysics.DENSITY,1)-node.GetSolutionStepValue(KratosMultiphysics.DENSITY))**2
                n_density     = n_density + (node.GetSolutionStepValue(KratosMultiphysics.DENSITY,1))**2   
                e_momentum  = e_momentum + (node.GetSolutionStepValue(KratosMultiphysics.MOMENTUM,1).norm_2()-node.GetSolutionStepValue(KratosMultiphysics.MOMENTUM).norm_2())**2
                n_momentum  = n_momentum + (node.GetSolutionStepValue(KratosMultiphysics.MOMENTUM,1).norm_2())**2 
                e_totalEnergy = e_totalEnergy + (node.GetSolutionStepValue(KratosMultiphysics.TOTAL_ENERGY,1)-node.GetSolutionStepValue(KratosMultiphysics.TOTAL_ENERGY))**2
                n_totalEnergy = n_totalEnergy + (node.GetSolutionStepValue(KratosMultiphysics.TOTAL_ENERGY,1))**2 
                
            e_density     = (e_density/n_density)**0.5
            e_momentum    = (e_momentum/n_momentum)**0.5
            e_totalEnergy = (e_totalEnergy/n_totalEnergy)**0.5
            
            if self.print_error_estimation:
                KratosMultiphysics.Logger.Print("               Density error estimation     : ", "{:.4e}".format(e_density)) 
                KratosMultiphysics.Logger.Print("               Momentum error estimation    : ", "{:.4e}".format(e_momentum)) 
                KratosMultiphysics.Logger.Print("               Total Energy error estimation: ", "{:.4e}".format(e_totalEnergy))

            if self.print_warnings:
                if math.isnan(e_totalEnergy):
                    KratosMultiphysics.Logger.Print("               NaN value in TOTAL ENERGY EQUATION")
                    sys.exit('               STOP')
                if math.isinf(e_totalEnergy):
                    KratosMultiphysics.Logger.Print("               Infinite value in TOTAL ENERGY EQUATION")
                    sys.exit('               STOP')
                if e_totalEnergy > 1e10:
                    KratosMultiphysics.Logger.Print("               Overflow in TOTAL ENERGY EQUATION")
                    KratosMultiphysics.Logger.Print("               ", e_totalEnergy)
                    sys.exit('               STOP')
                        
                if math.isnan(e_momentum):
                    KratosMultiphysics.Logger.Print("               NaN value in MOMENTUM X EQUATION")
                    sys.exit('               STOP')
                if math.isinf(e_momentum):
                    KratosMultiphysics.Logger.Print("               Infinite value in MOMENTUM X EQUATION")
                    sys.exit('               STOP')
                if e_momentum > 1e10:
                    KratosMultiphysics.Logger.Print("               Overflow in MOMENTUM EQUATION")
                    KratosMultiphysics.Logger.Print("               ", e_momentum)
                    sys.exit('               STOP')
                        
                if math.isnan(e_density):
                    KratosMultiphysics.Logger.Print("               NaN value in DENSITY EQUATION")
                    sys.exit('               STOP')
                if math.isinf(e_density):
                    KratosMultiphysics.Logger.Print("               Infinite value in DENSITY EQUATION")
                    sys.exit('               STOP')
                if e_density > 1e10:
                    KratosMultiphysics.Logger.Print("               Overflow in DENSITY EQUATION")
                    KratosMultiphysics.Logger.Print("               ", e_density)
                    sys.exit('               STOP')

            total_error = (e_density+e_totalEnergy+e_momentum)/3.0
            if total_error <= self.tolerance:
                KratosMultiphysics.Logger.Print("               Final total error:")
                KratosMultiphysics.Logger.Print("               ", total_error)
                sys.exit('               FINISHED')

            KratosMultiphysics.Logger.Print("::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
