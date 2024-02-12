#makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

import KratosMultiphysics, inspect, sys, timeit

from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis import PotentialFlowAnalysis
import KratosMultiphysics.CompressiblePotentialFlowApplication as KratosCompressFlow
import numpy as np

def inspect_class(myclass):
	'''Aux function to inspect all the methods in a class'''
	f=open("attributes.dat",'w')
	val=inspect.getmembers(myclass, lambda a:not(inspect.isroutine(a)))
	f.write("%s" %(val))
	f.close()

class PotentialFlowAnalysiswithOUTPUT(PotentialFlowAnalysis):
    '''Adaptation of PotentialFlowAnalysis class to output lift and cp'''

    def __init__(self,model,project_parameters):
        super(PotentialFlowAnalysiswithOUTPUT,self).__init__(model,project_parameters)

    def Finalize(self):
        super(PotentialFlowAnalysiswithOUTPUT,self).Finalize()
        # Printing cp
        self.write_cp("FluidModelPart.walls","walls.dat")
        # Printing lift
        self.write_lift("FluidModelPart.walls","coeffs.dat")

    def write_cp(self,mypartname,outpath):    
        # NOTE: Checking what attributes we have access to
        # inspect_class(KratosMultiphysics)
        mpart = self.model[mypartname]
        fout=open(outpath,'w')
        fout.write("#  CoordinatesX CoordinatesY CoordinatesZ CoefPressure\n")
        for node in mpart.Nodes:
            # NOTE: Checking what attributes we have access to
            #inspect_class(node)
            x=node.X0 ; y=node.Y0 ; z=node.Z0
            val=node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
            fout.write("%s %s %s %s\n" %(x,y,z,val))
        fout.close()
    
    def write_lift(self,mypartname,outpath): 
        mpart = self.model[mypartname]
        cl = mpart.ProcessInfo.GetValue(KratosCompressFlow.LIFT_COEFFICIENT)
        fout=open(outpath,'w')
        fout.write("%s" %(cl))
        fout.close()

if __name__ == "__main__":
    # We declare a user-defined time tracker
    wclockstart= timeit.default_timer()
    with open("ProjectParameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = PotentialFlowAnalysiswithOUTPUT(model,parameters)

    simulation.Run()

    # We write the total wall clock time to a file
    wclockend= timeit.default_timer()
    deltatime=wclockend-wclockstart
    fout=open("times.dat","w")
    fout.write("#wall clock time \n %s" %(deltatime))
    fout.close()
