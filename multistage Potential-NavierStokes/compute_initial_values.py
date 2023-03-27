# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as KratosCompressible


def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ComputeInitialValues(model, settings["Parameters"])

# All the processes python should be derived from "Process"
class ComputeInitialValues(KratosMultiphysics.Process):

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
                        "model_part"  : ""
        }""")

        settings.ValidateAndAssignDefaults(default_settings)
        self.model_part_name = settings["model_part"].GetString()
        # saving the modelpart
        self.model_part = model[self.model_part_name]
        print("::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 1")
        
    def ExecuteFinalizeSolutionStep(self):
        step = self.model_part.ProcessInfo[KratosMultiphysics.STEP]
        print("::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 2")
        if step == 1:
            nodal_variables = ["VELOCITY","PRESSURE"]
            pot_model_part = self.model_part
            nodal_value_process = KratosCompressible.ComputeNodalValueProcess(pot_model_part, nodal_variables)
            nodal_value_process.Execute()
            # Transfer the potential flow values as initial condition for the compressible problem
            tem_0 = 273
            gamma = 1.4
            c_v = 722.14
            c = (gamma*(gamma-1.0)*c_v*tem_0)**0.5
            free_stream_rho = 1.0
            free_stream_mach = 0.8
            pot_nodes = pot_model_part.Nodes
            comp_nodes = self._GetSolver().GetComputingModelPart().Nodes
            for pot_node, comp_node in zip(pot_nodes, comp_nodes):
                pot_v = pot_node.GetValue(KratosMultiphysics.VELOCITY)
                pot_v_norm_2 = pot_v[0] * pot_v[0] + pot_v[1] * pot_v[1]
                pot_v_norm = pot_v_norm_2**0.5
                mach = pot_v_norm / c
                num = 1.0 + 0.5 * (gamma - 1.0) * free_stream_mach**2
                det = 1.0 + 0.5 * (gamma - 1.0) * mach**2
                pot_rho = free_stream_rho * (num / det)**(1.0 / (gamma - 1.0))
                aux_tot_ener = pot_rho * (c_v * tem_0 + 0.5 * pot_v_norm_2)

                # Density initial condition
                comp_node.SetSolutionStepValue(KratosMultiphysics.DENSITY, 0, pot_rho)
                # Momentum initial condition
                comp_node.SetSolutionStepValue(KratosMultiphysics.MOMENTUM, 0, pot_rho * pot_v)
                # Total energy initial condition
                comp_node.SetSolutionStepValue(KratosMultiphysics.TOTAL_ENERGY, 0, aux_tot_ener)
