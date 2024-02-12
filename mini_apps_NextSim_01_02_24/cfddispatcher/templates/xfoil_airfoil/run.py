#!/usr/bin/env python
import timeit, sys
sys.path.append("$cfddfold")
import cfddispatcher as cfdd

########### INPUT PARAMETERS ############

# Mach number. If set to None, incompressible flow is assumed
M=$Mach

# Reynolds number. If set to None, inviscid computation is assumed
Re=$Re

# Angle of attack (degrees)
alpha=$alphadeg

# n_crit. If None, default value of 9 will be used. If set to 1, the run should correspond to a bypass transition
# NOTE: Only used for viscous runs
n_crit=$n_crit

# Values where transition will occur (top,bottom) if not trigered by n_crit 
#   If set to None, xfoil sets the default values to be [1.0,1.0])
# NOTE: Only used for viscous runs
xtr=$xtr

# Maximum number of iterations. Let it be None to keep the default of xfoil.
# NOTE: For initial inviscid tests, there is very low sensitivity with regards to this one.
max_iter=$max_iter

# File containing the coordinates of the airfoil
# NOTE: Should be provided starting from TE, and going either clockwise or counter clockwise
# NOTE: When dealing with close geometries (e.g. airfoils with sharp TEs), 
#       the TE point should be included twice
shapecoord="$shapecoord"

# Number of panels to use for re-paneling (required)
nrepanel=$nrepanel

# Deactivate the graphical option
# IMPORTANT: Set it to True. If you get an error message like "Program received signal SIGFPE: Floating-point exception - erroneous arithmetic operation." when
#			 trying to run xfoil, fall back into False 
nograph=$nograph

########### MAIN CODE ##################

# We declare a user-defined time tracker (based on wall clock time, so that
# other processes could interfere)
wclockstart= timeit.default_timer()
# 
# We perform our xfoil run
myx=cfdd.solvers_api.xfoilshape(shapecoord,nrepanel)
myx.set_max_iter(max_iter)
myx.set_n_crit(n_crit)
myx.set_xtr(xtr)
myx.set_nograph(nograph)
# NOTE: One could give the flag deletein=False to inspect the input
# NOTE: One could give the flag deleteout=False to keep all the standard
#       output files from xfoil. Not all the information is currently
#       stored in the class, so one would like to extend the scanning
#       functions for particular applications.
myx.run_xfoil(alpha, M, Re)
myx.get_coeffs()
myx.get_distrib()
myx.clean_output()

# We compute the total wall clock time
wclockend= timeit.default_timer()
deltatime=wclockend-wclockstart

# ........ Write files ............
fout=open("times.dat","w")
fout.write("#wall clock time[s] \n %s" %(deltatime))
fout.close()
#
fout=open("coeffs.dat",'w')
fout.write("#cl cd cm\n")
fout.write("%s %s %s" %(myx.CL,myx.CD,myx.CM))
fout.close()
#
fout=open("walls.dat",'w')
fout.write("# x y z cp\n")
for xnow,ynow,cpnow in zip(myx.xout,myx.yout,myx.Cp):
    fout.write("%s %s 0.0 %s\n" %(xnow,ynow,cpnow))
fout.close()
