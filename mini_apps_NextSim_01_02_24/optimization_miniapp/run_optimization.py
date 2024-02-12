import numpy as np
from shutil import rmtree
from geneticalgorithm import geneticalgorithm as ga
import os
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

### INTERACTION WITH MODELS ####
class fidelity:
  '''Allows to switch the use of lower or higher fidelity models in the numerical optimization'''
  def __init__(self,pyvers,scriptname,local_fold,plabel,m,p,Mach,adeg): 
    self.pyvers=pyvers
    self.scriptname=scriptname
    self.local_fold=local_fold
    self.plabel=plabel
    #
    self.m=m ; self.p=p
    self.Mach=Mach ; self.adeg=adeg
    #
    self.simfold=os.path.join(local_fold,plabel)
    self.histvar=[] # storing design variables for every call
    self.histobj=[] # storing objective functions for every call
  def run_model(self,t):
    '''Takes as an input the design variable (thickess), runs a simulation, 
    and returns the objective function'''
    self.run_simulation(t)
    self.get_cl()
    self.clean_simulation()
    self.histvar.append(t)
    self.histobj.append(self.cl)
    return -self.cl
  def run_simulation(self,t):
    '''Runs the right model'''
    runcmd="python %s.py %s %s %s %s %s %s %s %s" %(self.scriptname,self.pyvers,
                                                  self.local_fold,self.plabel,t[0],
                                                  self.m,self.p,self.Mach,self.adeg)
    os.system(runcmd)
  def get_cl(self):
    coeffpath=os.path.join(self.simfold,"coeffs.dat")
    self.cl=np.loadtxt(coeffpath,usecols=(0,))
  def clean_simulation(self):
    rmtree(self.simfold)

######## USER-RELATED INPUT DATA ########
# Python version to execute optimization
pyvers="python3.10"

###### INPUT DATA ######
# Airfoil geometry
m=7.0  # m/100, max chamber.
p=4.0  # p/10.0: position of max chamber
# Inflow
Mach=0.5 # Mach number
adeg=3.0  # Angle of attack in degrees
varbound=np.array([[10,20]]) # Definition of airfoil thickness (/100.0) bounds

######## OPTIMIZATION SETTINGS ####
maxoptiter=8 # Max. number of optimization iterations
popsize=8 # Population size
maxnoimprov=maxoptiter # Maximum iterations without improvement

######## SCRIPT PARAMETERS #####
# Folder to store results
local_fold="results"
# label for our simulations
plabel="sim"

####### INSTANTIATING MODELS FOR OPTIMIZATION ######
fidHF=fidelity(pyvers,"full_potential",local_fold,plabel,
               m,p,Mach,adeg) # Higher fidelity model 
fidLF=fidelity(pyvers,"panel",local_fold,plabel,
               m,p,Mach,adeg) # Lower fidelity model 

####### PERFORMING OPTIMIZATIONS ######
opt_dict={'max_num_iteration': maxoptiter,
          'population_size':popsize,
          'mutation_probability':0.1,
          'elit_ratio': 0.01,
          'crossover_probability': 0.5,
          'parents_portion': 0.3,
          'crossover_type':'uniform',
          'max_iteration_without_improv':maxnoimprov}
optLF=ga(function=fidLF.run_model,dimension=1,variable_type='real',variable_boundaries=varbound,
         function_timeout=5*60,
         algorithm_parameters=opt_dict,
         convergence_curve=False,progress_bar=False)
optHF=ga(function=fidHF.run_model,dimension=1,variable_type='real',variable_boundaries=varbound,
         function_timeout=5*60,
         algorithm_parameters=opt_dict,
         convergence_curve=False,progress_bar=False)
optLF.run()
optHF.run()

############ PLOTTING ################## 
fig, ax = plt.subplots()
fig.set_figwidth(9.0)
fig.set_figheight(3.4)
#
ax.plot(fidLF.histvar,fidLF.histobj,
  "bo",markersize=6.0,label="Simulations LF")
ax.plot(optLF.output_dict["variable"],-optLF.output_dict["function"],
  "rx",markersize=6.0,label="Optimum LF")
ax.plot(fidHF.histvar,fidHF.histobj,
  "gs",markersize=6.0,label="Simulations HF")
ax.plot(optHF.output_dict["variable"],-optHF.output_dict["function"],
  "r+",markersize=6.0,label="Optimum HF")
#
ax.grid()
plt.legend()
plt.tight_layout()
plt.xlabel("t")
plt.ylabel("cl")
fig.savefig("optima_found.png")
plt.close("all")
