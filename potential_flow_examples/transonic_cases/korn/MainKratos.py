import sys
import time
import importlib
import numpy as np
import matplotlib.pyplot as plt
# import matplotlib
# matplotlib.use('Agg')

import KratosMultiphysics

def CreateAnalysisStageWithFlushInstance(cls, global_model, parameters):
    class AnalysisStageWithFlush(cls):

        def __init__(self, model,project_parameters, flush_frequency=10.0):
            super().__init__(model,project_parameters)
            self.flush_frequency = flush_frequency
            self.last_flush = time.time()
            sys.stdout.flush()

        def Initialize(self):
            super().Initialize()
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

    with open("ProjectParameters.json", 'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    analysis_stage_module_name = parameters["analysis_stage"].GetString()
    analysis_stage_class_name = analysis_stage_module_name.split('.')[-1]
    analysis_stage_class_name = ''.join(x.title() for x in analysis_stage_class_name.split('_'))

    analysis_stage_module = importlib.import_module(analysis_stage_module_name)
    analysis_stage_class = getattr(analysis_stage_module, analysis_stage_class_name)
    
    global_model = KratosMultiphysics.Model()
    simulation = CreateAnalysisStageWithFlushInstance(analysis_stage_class, global_model, parameters)
    simulation.Run()

    modelpart = global_model["FluidModelPart.Body2D_Body"]
    x = np.zeros(modelpart.NumberOfNodes())
    y = np.zeros(modelpart.NumberOfNodes())
    z = np.zeros(modelpart.NumberOfNodes())
    cp = np.zeros(modelpart.NumberOfNodes())
    for i,node in enumerate(modelpart.Nodes):
        x[i] = node.X0 ; y[i] = node.Y0 ; z[i] = node.Z0
        cp[i] = node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
    
    # Plot cp vs x
    fig,ax  = plt.subplots()
    fig.set_figwidth(10.0)
    fig.set_figheight(6.0)
    ax.plot( x, cp, "o", markersize = 5.0)
    ax.grid()
    plt.ylabel('Cp')
    plt.xlabel('x')
    plt.title('Korn’s supercritical airfoil -  Alpha = 0.7º Mach = 0.75')
    ax.invert_yaxis()
    plt.tight_layout()
    fig.savefig("Airfoils_Cp_x.png", dpi=400)
    plt.close('all')

    
    # from matplotlib import cbook
    # fig, ax = plt.subplots()
    # fig.set_figwidth(10.0)
    # fig.set_figheight(6.0)
    # # make data
    # ax.plot( x, -cp, "o", markersize = 5.0)
    # # ax.invert_yaxis()
    # ax.grid()
    # # fig = plt.axis([-0.55,1.5,np.max(cp)+0.1,np.min(cp)-0.1])
    # fig = plt.ylabel('Cp')
    # fig = plt.xlabel('x')
    # fig = plt.title('Korn’s supercritical airfoil -  Alpha = 0.7º Mach = 0.75')
    # # inset Axes....
    # x1, x2, y1, y2 = 0.49, 0.51, -0.5, -0.2  # subregion of the original image
    # axins = ax.inset_axes( [0.4, 0.15, 0.25, 0.25],
    #     xlim=(x1, x2), ylim=(y1, y2), xticklabels=[], yticklabels=[])
    # axins.plot( x, -cp, ".-")
    # # axins.invert_yaxis()
    # axins.grid()
    # ax.indicate_inset_zoom(axins, edgecolor="black")
    # fig = plt.tight_layout()
    # fig = plt.savefig("Airfoils_Cp_x_zoom.png", dpi=400)
    # fig = plt.close('all')
