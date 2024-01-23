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
                                "print_max_values"       : false,
                                "variable_max_value"     : "",
                                "print_error_estimation" : false,
                                "print_warnings"         : false,
                                "model_part"             : "FluidModelPart.fluid_computational_model_part"
        }""")

        settings.ValidateAndAssignDefaults(default_settings)
        self.print_max_values = settings["print_max_values"].GetBool()
        self.variable_max_value = settings["variable_max_value"].GetString()
        self.print_warnings = settings["print_warnings"].GetBool()
        self.print_error_estimation = settings["print_error_estimation"].GetBool()
        self.model_part_name = settings["model_part"].GetString()
                
        # saving the modelpart
        self.model_part = model[self.model_part_name]
    
    def ExecuteBeforeSolutionLoop(self):
        if self.print_max_values or self.print_warnings or self.print_error_estimation:
            KratosMultiphysics.Logger.Print("::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
            KratosMultiphysics.Logger.Print("::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
            
    def ExecuteFinalizeSolutionStep(self):
        if self.print_max_values or self.print_warnings or self.print_error_estimation:
            if self.print_max_values:
                max_density = 0.0
                max_energy = 0.0
                max_momentum_module = 0.0
                for node in self.model_part.Nodes:
                    density  = node.GetSolutionStepValue(KratosMultiphysics.DENSITY)
                    energy   = node.GetSolutionStepValue(KratosMultiphysics.TOTAL_ENERGY)
                    momentum = node.GetSolutionStepValue(KratosMultiphysics.MOMENTUM).norm_2()
                    if max_density < density :
                        max_density = density
                    if max_energy < energy :
                        max_energy = energy
                    if max_momentum_module < momentum :
                        max_momentum_module = momentum

                if self.variable_max_value == "density":
                    KratosMultiphysics.Logger.Print("               Max Density  : ", max_density)
                elif self.variable_max_value == "energy":
                    KratosMultiphysics.Logger.Print("               Max Energy   : ", max_energy)
                elif self.variable_max_value == "momentum":
                    KratosMultiphysics.Logger.Print("               Max Momentum : ", max_momentum_module)
                else:
                    KratosMultiphysics.Logger.Print("               Max Density  : ", max_density)
                    KratosMultiphysics.Logger.Print("               Max Energy   : ", max_energy)
                    KratosMultiphysics.Logger.Print("               Max Momentum : ", max_momentum_module)
                
            step = self.model_part.ProcessInfo.GetValue(KratosMultiphysics.STEP)
            if  step > 1 :
                e_density     = 0  
                e_momentum_x  = 0
                e_momentum_y  = 0
                e_totalEnergy = 0
                n_density     = 0  
                n_momentum_x  = 0
                n_momentum_y  = 0
                n_totalEnergy = 0
                for node in self.model_part.Nodes:
                    e_density     = e_density + (node.GetSolutionStepValue(KratosMultiphysics.DENSITY,1)-node.GetSolutionStepValue(KratosMultiphysics.DENSITY))**2
                    n_density     = n_density + (node.GetSolutionStepValue(KratosMultiphysics.DENSITY,1))**2   
                    e_momentum_x  = e_momentum_x + (node.GetSolutionStepValue(KratosMultiphysics.MOMENTUM_X,1)-node.GetSolutionStepValue(KratosMultiphysics.MOMENTUM_X))**2
                    n_momentum_x  = n_momentum_x + (node.GetSolutionStepValue(KratosMultiphysics.MOMENTUM_X,1))**2 
                    e_momentum_y  = e_momentum_y + (node.GetSolutionStepValue(KratosMultiphysics.MOMENTUM_Y,1)-node.GetSolutionStepValue(KratosMultiphysics.MOMENTUM_Y))**2
                    n_momentum_y  = n_momentum_y + (node.GetSolutionStepValue(KratosMultiphysics.MOMENTUM_Y,1))**2 
                    e_totalEnergy = e_totalEnergy + (node.GetSolutionStepValue(KratosMultiphysics.TOTAL_ENERGY,1)-node.GetSolutionStepValue(KratosMultiphysics.TOTAL_ENERGY))**2
                    n_totalEnergy = n_totalEnergy + (node.GetSolutionStepValue(KratosMultiphysics.TOTAL_ENERGY,1))**2 
                   
                e_density     = (e_density/n_density)**0.5
                e_momentum_x  = (e_momentum_x/n_momentum_x)**0.5
                e_momentum_y  = (e_momentum_y/n_momentum_y)**0.5
                e_totalEnergy = (e_totalEnergy/n_totalEnergy)**0.5
                    
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
                            
                    if math.isnan(e_momentum_x):
                        KratosMultiphysics.Logger.Print("               NaN value in MOMENTUM X EQUATION")
                        sys.exit('               STOP')
                    if math.isinf(e_momentum_x):
                        KratosMultiphysics.Logger.Print("               Infinite value in MOMENTUM X EQUATION")
                        sys.exit('               STOP')
                    if e_momentum_x > 1e10:
                        KratosMultiphysics.Logger.Print("               Overflow in MOMENTUM X EQUATION")
                        KratosMultiphysics.Logger.Print("               ", e_momentum_x)
                        sys.exit('               STOP')
                          
                    if math.isnan(e_momentum_y):
                        KratosMultiphysics.Logger.Print("               NaN value in MOMENTUM Y EQUATION")
                        sys.exit('               STOP')
                    if math.isinf(e_momentum_y):
                        KratosMultiphysics.Logger.Print("               Infinite value in MOMENTUM Y EQUATION")
                        sys.exit('               STOP')
                    if e_momentum_y > 1e10:
                        KratosMultiphysics.Logger.Print("               Overflow in MOMENTUM Y EQUATION")
                        KratosMultiphysics.Logger.Print("               ", e_momentum_y)
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
                    
                if self.print_error_estimation:
                    KratosMultiphysics.Logger.Print("               Density error estimation     : ", "{:.4e}".format((e_density/n_density)**0.5)) 
                    KratosMultiphysics.Logger.Print("               Momentum_x error estimation  : ", "{:.4e}".format((e_momentum_x/n_momentum_x)**0.5)) 
                    KratosMultiphysics.Logger.Print("               Momentum_Y error estimation  : ", "{:.4e}".format((e_momentum_y/n_momentum_y)**0.5)) 
                    KratosMultiphysics.Logger.Print("               Total Energy error estimation: ", "{:.4e}".format((e_totalEnergy/n_totalEnergy)**0.5))
                    
            KratosMultiphysics.Logger.Print("::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
            KratosMultiphysics.Logger.Print("::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
