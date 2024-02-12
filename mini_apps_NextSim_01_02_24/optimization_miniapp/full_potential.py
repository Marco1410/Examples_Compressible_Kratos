import sys,os
sys.path.append('../')
import cfddispatcher as cfdd
import numpy as np

####################################
# USER-RELATED INPUT DATA ##########
####################################
# Parameters for local execution 
mykratosfold="/home/sgho/CIMNE_Documents/11_FULL_POTENTIAL_KRATOS/021_Kratos_installation_persoHP_laptop"
kratosenvdict={
"PYTHONPATH": "%s/Kratos/bin/Release" %(mykratosfold),
"LD_LIBRARY_PATH": "%s/Kratos/bin/Release/libs" %(mykratosfold)
}
# Limit of threads to use when running Kratos
# NOTE: Let it be None to avoid imposing any limit
kratosnthread_lim=12

##############################
## INPUT VARIABLES ##########
##############################
kratosexec=str(sys.argv[1]) # Python version to execute simulation
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
nairfoilsurf=500
# Size of CFD domain
widthL=20.0
widthR=70.0
heightD=25.0
heightU=25.0
# Mesh paramaters
extsize=2.0 # outer domain
lcwall=0.00625 # wall refinement
lcedgemin=0.00375 # wall refinement near LE
lcedgemax=0.00625 # wall refinement near TE
distmax=20.0 # distance from wall to relax to extsize
blfirst=lcwall # first cell height for extruded BL
blthick=8*lcwall # total height of extruded BL
blratio=1.0 # expansion ration of extruded BL
nfanmax=25 # number of fans at TE
# Kratos model to use
tfolderkratos="../cfddispatcher/templates/kratos_airfoil_compressible"
# Kratos max iterations per set of params
kratosmaxel=3000

#####################################
# RUN SIMULATION 
#####################################
# .................
# Initialize
# .................
localsimgroup=cfdd.links.localsimgroup(locfold=local_fold)
# .................
# Generate the mesh
# .................
myparam=cfdd.param.NACAparam(m=m,p=p,t=t)
mygeom=cfdd.geom.airfoil(param=myparam,TEtype="sharp")
mygeom.generate_geometry(n=nairfoilsurf,distri="sine", redist=True, rechord=True)
#
kratosmeshbase="mesh_%s" %(plabel)
kratosmeshpath="%s.msh" %(kratosmeshbase)
kratosmmdpapath="%s.mdpa" %(kratosmeshbase)
xl,yl=mygeom.give_coordinates("PS")
xu,yu=mygeom.give_coordinates("SS")
mymesh=cfdd.mesh.extaero(widthL=widthL,widthR=widthR,
                         heightD=heightD,heightU=heightU,
                         xl=xl,yl=yl,xu=xu,yu=yu)
mymesh.generate_geom_and_mesh(kratosmeshpath,
						  lcedgemin=lcedgemin,lcedgemax=lcedgemax,
						  lcwall=lcwall,extsize=extsize,
						  blthick=blthick,blfirst=blfirst,blratio=blratio,  
						  nfanmax=nfanmax,  
                          threshold={"distmin":blthick,"distmax":distmax,"lcmin":lcwall},
						  topology="tri",
                          planeid="xy")
cfdd.mesh.from_gmsh41_to_Kratosmdpa(kratosmeshpath,kratosmmdpapath)
# .......................
# Add solver simulation
# .......................
tdict={}
tdict["rho"]="1.225"
tdict["heatratio"]="1.4"
tdict["MaxNumIterations"]="%s" %(kratosmaxel)
tdict["angleOfAttack"]="%s" %(np.deg2rad(adeg))
tdict["Mach"]="%s" %(Mach)
tdict["speed_of_sound"]="340.0"
tdict["meshbasename"]=kratosmeshbase
tdict["simlabel"]=plabel
tdict["removemesh"]="True" # Clean up mesh after finishing Kratos computation
# 	
mysim=cfdd.links.sim(tfolderkratos,tdict)
mysim.add_input_data(kratosmmdpapath)
#
launcher=cfdd.queues.mypc()
run_commandkratos='%s MainKratos.py > konsole.out' %(kratosexec)
launcher.define_run_command(run_commandkratos)
mysim.add_launcher(launcher)
#
localsimgroup.add_simulation(mysim,plabel)
# .......................
# Run simulation
# .......................
localsimgroup.deploy_sims_locally()
localsimgroup.copy_additional_data(keepsourcedata=False)
if(kratosnthread_lim!=None):
	kratosenvdict["OMP_NUM_THREADS"]=kratosnthread_lim
localsimgroup.perform_simulations(kratosenvdict)

