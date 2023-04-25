import KratosMultiphysics

class PrintModelModeler(KratosMultiphysics.Modeler):

    def __init__(self, model, settings):
        super().__init__(model, settings)

        # Cannot validate as settings may differ among input types
        settings.AddMissingParameters(self.__GetDefaultSettings())

        # Declare required member variables
        self.model = model
        self.settings = settings

    def SetupGeometryModel(self):
        super().SetupGeometryModel()

    def PrepareGeometryModel(self):
        super().PrepareGeometryModel()

    def SetupModelPart(self):
        super().SetupModelPart()

        output_filename = self.settings["output_filename"].GetString()
        model_part_name = self.settings["model_part_name"].GetString()
        model_part = self.model.GetModelPart(model_part_name)
        KratosMultiphysics.ModelPartIO(output_filename, KratosMultiphysics.IO.WRITE | KratosMultiphysics.IO.MESH_ONLY).WriteModelPart(model_part)

    def __GetDefaultSettings(self):
        default_settings = KratosMultiphysics.Parameters('''{
                "output_filename" : "",
                "model_part_name" : ""
        }''')
        return default_settings

def Factory(model, settings):
    return PrintModelModeler(model, settings)

