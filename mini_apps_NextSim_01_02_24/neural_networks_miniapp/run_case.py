import sys,os
sys.path.append('../')
import cfddispatcher as cfdd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

print("+++++++++++ START RUNNING run_case ++++++++++++++++")

####################################
### INPUT DATA #####################
####################################
# Airfoil geometry
m=3.0  # m/100, max chamber. Training done for [2.0,5.0]
p=4.0  # p/10.0: position of max chamber. Training done for [3.0,6.0]
t=11.5 # t/100.0: thickness. Training done for [11.0,12.0]
# Inflow
Mach=0.68  # Mach number. Training done for [0.65,0.7]
adeg=0.5   # Angle of attack in degrees. Training done for [0.0,1.0]

####################################
# USER-RELATED INPUT DATA ##########
####################################
# Parameters for local execution 
kratosexec="python3.10"
mykratosfold="../../"
kratosenvdict={
"PYTHONPATH": "%s/Kratos/bin/Release" %(mykratosfold),
"LD_LIBRARY_PATH": "%s/Kratos/bin/Release/libs" %(mykratosfold)
}
# Limit of threads to use when running Kratos
# NOTE: Let it be None to avoid imposing any limit
kratosnthread_lim=12

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
tfolderkratos="../cfddispatcher/templates/kratos_airfoil_transonic"
# Kratos max iterations per set of params
# NOTE: kept very low as we will use the outer loop with several iterations
kratosmaxel=50

###############################
## SCRIPT PARAMETERS ##########
###############################
# Folder to store results
local_fold="results"
# Label for simulation
plabel="KratosFpotential"

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

#####################################
# GET RESULTS 
#####################################
cpfold=os.path.join(local_fold,plabel)
cppath=os.path.join(cpfold,"walls.dat")
if(not os.path.isfile(cppath)):
	print("ERROR: something went wrong with Kratos simulation")
	exit()
data = np.loadtxt(cppath)
#####################################
# APPLY CORRECTION 
#####################################
datacorr=data.copy() 
# TODO: Just applying a factor here for illustration, but we should call the correction NN
datacorr[:,3]= datacorr[:,3]*1.1
#####################################
# PLOT COMPARISON
#####################################
fig, ax = plt.subplots()
fig.set_figwidth(9.0)
fig.set_figheight(3.4)
#
ax.plot(data[:,0],data[:,3],
	"ko",markersize=2.0,label="Full Potential (FP)")
ax.plot(datacorr[:,0],datacorr[:,3],
	"rx",markersize=2.0,label="FP + Euler correction")
#
ax.grid()
plt.legend()
plt.tight_layout()
plt.gca().invert_yaxis()
plt.xlabel("x")
plt.ylabel("cp")
fig.savefig("pressure_coeffs.png")
plt.close("all")

print("+++++++++++ END RUNNING run_case ++++++++++++++++")
