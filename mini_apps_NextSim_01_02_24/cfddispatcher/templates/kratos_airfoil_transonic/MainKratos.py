#makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

import KratosMultiphysics, inspect, sys, timeit, os

import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis import PotentialFlowAnalysis
import numpy as np

class PotentialFlowAnalysisWithFlush(PotentialFlowAnalysis):
    '''Adaptation of PotentialFlowAnalysis class to 1) make iterative runs and 2) output lift and cp'''

    def __init__(self,model,project_parameters):
        super(PotentialFlowAnalysisWithFlush,self).__init__(model,project_parameters)

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

    def ModifyAfterSolverInitialize(self):
        if (os.path.isfile("previous_data_1_2.dat")):
            id_nodes = np.loadtxt("previous_data_1_2.dat",usecols=(0,))
            velocity_potential = np.loadtxt("previous_data_1_2.dat",usecols=(1,))
            auxiliary_velocity_potential = np.loadtxt("previous_data_1_2.dat",usecols=(2,))
            for i in range(len(id_nodes)):
                node = model["FluidModelPart"].GetNode(int(id_nodes[i]))
                node.SetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL, 0, velocity_potential[i])
                node.SetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL, 1, velocity_potential[i])
                node.SetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL, 0, auxiliary_velocity_potential[i])
                node.SetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL, 1, auxiliary_velocity_potential[i])
        elif(os.path.isfile("previous_data_1.dat")):
            id_nodes = np.loadtxt("previous_data_1.dat",usecols=(0,))
            velocity_potential = np.loadtxt("previous_data_1.dat",usecols=(1,))
            auxiliary_velocity_potential = np.loadtxt("previous_data_1.dat",usecols=(2,))
            for i in range(len(id_nodes)):
                node = model["FluidModelPart"].GetNode(int(id_nodes[i]))
                node.SetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL, 0, velocity_potential[i])
                node.SetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL, 1, velocity_potential[i])
                node.SetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL, 0, auxiliary_velocity_potential[i])
                node.SetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL, 1, auxiliary_velocity_potential[i])

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

    def Finalize(self):
        super(PotentialFlowAnalysisWithFlush,self).Finalize()
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
        cl = mpart.ProcessInfo.GetValue(CPFApp.LIFT_COEFFICIENT)
        fout=open(outpath,'w')
        fout.write("%s" %(cl))
        fout.close()

def M_cr(cp,critical_mach):
    return cp-((1/(0.7*critical_mach**2))*((((2+0.4*critical_mach**2)/2.4)**3.5)-1))

def transonic_parameters(mach_infinity):
    cp_min = -0.6
    if mach_infinity < 0.7:
        # Prandtl-Glauert rule (not precise for M > 0.7)
        cp_M = cp_min/(np.sqrt(1-mach_infinity**2))
    else:
        # Karman-Tsien rule
        cp_M = cp_min/(np.sqrt(1-mach_infinity**2)+(cp_min*mach_infinity**2/(1+np.sqrt(1-mach_infinity**2))))
    # Bisection method
    a = 0.02
    b = 0.95
    tol=1.0e-6
    critical_mach = (a + b) / 2.0
    while True:
        if b - a < tol:
            return critical_mach
        elif M_cr(cp_M,a) * M_cr(cp_M,critical_mach) > 0:
            a = critical_mach
        else:
            b = critical_mach
        critical_mach = (a + b) / 2.0

def upwind_aprox(angle,mach):
    '''Angle of attack in degrees'''
    # interpolated values, using NACA0012, in the range: 
    # angle: 0.0, 5.0
    # mach : 0.60, 0.85
    nacaupw=2.1489732 /(1 + np.exp(-(1.55470776*angle+37.876707*mach-31.48703546))) + 0.59767296
    # IMPORTANT: If solution is very unstable, switch to a fix value of 0.1 
    # nacaupw=0.1
    return nacaupw

def inspect_class(myclass):
	'''Aux function to inspect all the methods in a class'''
	f=open("attributes.dat",'w')
	val=inspect.getmembers(myclass, lambda a:not(inspect.isroutine(a)))
	f.write("%s" %(val))
	f.close()

def remove_results(label):
    for filein in ["walls.dat","coeffs.dat"]:
        if(os.path.isfile(filein)):
            print("INFO: deleted %s file as %s" %(filein,label))
            os.remove(filein) 

if __name__ == "__main__":
    # We declare a user-defined time tracker (based on wall clock time, so that
    # other processes could interfere)
    wclockstart= timeit.default_timer()

    with open("ProjectParameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    # Prepare for the iterative procedure
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('previous_data_1.dat')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('previous_data_1_2.dat')
 
    mach_infinity=$Mach
    angledeg=np.rad2deg($angleOfAttack)
    critical_mach=transonic_parameters(mach_infinity)
    upwind_factor_constant = upwind_aprox(angledeg, mach_infinity)
    #
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach_infinity)
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["critical_mach"].SetDouble(critical_mach)
    #
    convergence_ratio = 1.0
    absolute_norm = 1.0
    tolerancia = 1e-10
    write_and_plot = True
    #
    while (convergence_ratio > tolerancia and absolute_norm > tolerancia):
        parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["upwind_factor_constant"].SetDouble(upwind_factor_constant)
        # Make sure last result files are deleted
        remove_results("INTERMEDIATE RESULTS")
        model = KratosMultiphysics.Model()
        simulation = PotentialFlowAnalysisWithFlush(model,parameters)
        simulation.Run()
        #
        convergence_ratio = model["FluidModelPart"].ProcessInfo[KratosMultiphysics.CONVERGENCE_RATIO]
        absolute_norm = model["FluidModelPart"].ProcessInfo[KratosMultiphysics.RESIDUAL_NORM]
        if convergence_ratio >= 0.5:
            upwind_factor_constant += 0.2
            #print("::::::::::::::::::::::::::::::::NEW UPWIND::::::::::::::::::::::::::::::::::::::")
        else:
            upwind_factor_constant += 0.05
            KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('previous_data_1_2.dat')
            fout=open("previous_data_1_2.dat",'w')
            modelpart = model["FluidModelPart"]
            for node in modelpart.Nodes:
                id_node = node.Id
                velocity_potential = node.GetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL,1)
                auxiliary_velocity_potential = node.GetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL,1)
                fout.write("%s %s %s\n" %(id_node, velocity_potential, auxiliary_velocity_potential))
            fout.close()
            #print(":::::::::::::::::::::::::::::::::RESTART / SAVE:::::::::::::::::::::::::::::::::")
        print("INTERATION KRATOS: tolerance %s, convergence_ratio %s, absolute_norm %s" %(tolerancia, convergence_ratio, absolute_norm))
        if upwind_factor_constant > 5.0:
            print("::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
            print("::::::::::::::::::::::::::::::::::::::: Non Converged ::::::::::::::::::::::::::::")
            print("::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
            write_and_plot = False
            # Make sure we do not provide results
            remove_results("NOT CONVERGED")
            break
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('previous_data_1_2.dat')
    if write_and_plot:
        fout=open("previous_data_1.dat",'w')
        modelpart = model["FluidModelPart"]
        for node in modelpart.Nodes:
            id_node = node.Id
            velocity_potential = node.GetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL,1)
            auxiliary_velocity_potential = node.GetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL,1)
            fout.write("%s %s %s\n" %(id_node, velocity_potential, auxiliary_velocity_potential))
        fout.close()

    # Clean up intermediate files related to the iterations
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('previous_data_1.dat')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('previous_data_1_2.dat')

    # Clean up mesh if requested by user
    if($removemesh):
        KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('$meshbasename.mdpa')

    # We write the total wall clock time to a file
    wclockend= timeit.default_timer()
    deltatime=wclockend-wclockstart
    fout=open("times.dat","w")
    fout.write("#wall clock time \n %s" %(deltatime))
    fout.close()
