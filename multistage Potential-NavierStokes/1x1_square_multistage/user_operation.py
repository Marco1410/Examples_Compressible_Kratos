import KratosMultiphysics as KMP
import KratosMultiphysics.python_operation as python_operation

class UserOperation(python_operation.PythonOperation):
    def __init__(self, settings, model):
        super().__init__()

    def Execute(self):
        KMP.Logger.PrintWarning("UserOperation", "Calling the fake UserOperation Execute().")

    @staticmethod
    def Create(settings, model):
        return UserOperation(settings, model)