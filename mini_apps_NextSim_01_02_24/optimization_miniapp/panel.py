import sys,os
sys.path.append('../')
import cfddispatcher as cfdd
import numpy as np

##############################
## INPUT VARIABLES ##########
##############################
xfoilexec=str(sys.argv[1]) # Python version to execute simulation
local_fold=str(sys.argv[2]) # folder to store results
plabel=str(sys.argv[3])  # label for simulation
t=float(sys.argv[4]) # airfoil thickness (DESIGN VARIABLE).
m=float(sys.argv[5]) # airfoil thickness (FIXED).
p=float(sys.argv[6]) # airfoil thickness (FIXED).
Mach=float(sys.argv[7]) # airfoil thickness (FIXED).
adeg=float(sys.argv[8]) # airfoil thickness (FIXED).

###############################
## MODEL PARAMETERS  ##########
###############################
# Number of points to generate airfoil
nairfoilsurf=200
# Xfoil model to use
tfolderxfoil="../cfddispatcher/templates/xfoil_airfoil"
# Where to store computed airfoil
shapecoord="shapecoord.csv"
#
#####################################
# RUN SIMULATION 
#####################################
# .................
# Initialize
# .................
localsimgroup=cfdd.links.localsimgroup(locfold=local_fold)
# ..............................
# Generate the airfoil points
# ...............................
myparam=cfdd.param.NACAparam(m=m,p=p,t=t)
mygeom=cfdd.geom.airfoil(param=myparam,TEtype="sharp")
mygeom.generate_geometry(n=nairfoilsurf,distri="sine", redist=True, rechord=True)
mygeom.write_airfoil(outfile=shapecoord,surftype="all")
# .......................
# Add solver simulation
# .......................
tdict={}
tdict["shapecoord"]=shapecoord
tdict["simlabel"]=plabel
tdict["Mach"]="%s" %(Mach)
tdict["Re"]="None"
tdict["alphadeg"]="%s" %(adeg)
tdict["n_crit"]="None"
tdict["xtr"]="None"
tdict["max_iter"]="500"
tdict["nrepanel"]="150"
# IMPORTANT: If you get an error message like "Program received signal SIGFPE: Floating-point exception - 
# erroneous arithmetic operation." when trying to run xfoil, fall back into False 
tdict["nograph"]="True"
#
mysim=cfdd.links.sim(tfolderxfoil,tdict)
mysim.add_input_data(shapecoord)
launcher=cfdd.queues.mypc()
run_command='%s run.py > konsole.out' %(xfoilexec)
launcher.define_run_command(run_command)
mysim.add_launcher(launcher)
localsimgroup.add_simulation(mysim,plabel)
# .......................
# Run simulation
# .......................
localsimgroup.deploy_sims_locally()
localsimgroup.copy_additional_data(keepsourcedata=False)
localsimgroup.perform_simulations()
