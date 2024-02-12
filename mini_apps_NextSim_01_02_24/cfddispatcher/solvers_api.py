#!/usr/bin/env python
import os,warnings,subprocess,time
import numpy as np

###################################################
# CLASSES WITH SOLVERS' API'S #####################
###################################################

class xfoilshape:
    ''' Class to prepare and perform xfoil simulations for single angles of attack.
    Results are stored in the class itself, after executing run_xfoil
    '''
    def __init__(self,infile,nrepanel):
        ''' Initializes xfoil computation based on geometry. 
            - infile contains teh 2d geometry (by columns, without delimiters, may have # comments)
            - nrepanel is an integer with the number of panels 
              that xfoil will create. Xfoil will always repanel in this class.
              Default set to Xfoil default.
        '''
        self.nrepanel=nrepanel
        if self.nrepanel==None:
	        warnings.warn('Remember that the number of points describing the airfoil should be smaller than 365 if repanelling is not used.')
        # ...... Variables initizialization ......
        self.max_iter=None
        self.n_crit=None
        self.xtr=None
        # ...... Parameters ......
        # Names of xfoil in files
        self.geom=infile
        self.rundef="rundef.in"
        # Names of xfoil out files
        self.coeffs="out_coeffs.out"
        self.cpfile="out_cps.out"
        self.extinfo="out_extended.out"
        # .... Clean in/out ...
        self._clean_input()
        self.clean_output()
    def set_max_iter(self,max_iter):
        ''' Sets the number of maximum iterations to perform '''
        self.max_iter=max_iter
    def set_n_crit(self,n_crit):
        ''' Sets n criticial for our simulations '''
        self.n_crit=n_crit
    def set_xtr(self,xtr):
        ''' Sets xtr for our simulations '''
        self.xtr=xtr
    def set_nograph(self,nograph):
        '''Sets the nograph flag, so that the graphical capabilities of xfoil can be deactivated'''
        self.nograph=nograph
    def run_xfoil(self, alphadeg, Mach=None, Re=None,
        deletein=True):
        ''' Runs computation for a given angle (in degrees),
            Mach number is optional (assumed incompressible flow if None). 
            Reynolds number is optional (assumed inviscid simulation if None)
            - deletin is a flag for removing xfoil input files
            - deleteout is a flag to remove xfoil output files 
        '''
        self.alphadeg=alphadeg
        self.Mach=Mach
        self.Re=Re
        # Create input in temporary files
        self._write_input_file()
        # run code 
        self._run_code()
        # Clean files if requested
        if(deletein==True):
            self._clean_input()
    def get_coeffs(self):
        ''' Gets aero coefficients once simulation is performed.
        Stores them in variables: CL,CD,CDp,CM,Top_Xtr,Bot_Xtr,Top_Itr,Bot_Itr''' 
        # Initialization of the variables
        self.CL=0.0 ; self.CD=0.0 ; self.CDp=0.0 ; self.CM=0.0
        self.Top_Xtr=0.0 ; self.Bot_Xtr=0.0 
        self.Top_Itr=0.0 ; self.Bot_Itr=0.0
        self.CL_CD=0.0
        #
        self.isconv=True
        if not os.path.isfile(self.coeffs):
            self.isconv=False
        else:
            with open(self.coeffs) as f:
                for line in f:
                    pass
                last_line = line
            sp=line.split()
            self.CL=sp[1] ; self.CD=sp[2] ; self.CDp=sp[3] ; self.CM=sp[4]
            self.Top_Xtr=sp[5] ; self.Bot_Xtr=sp[6] 
            if self.Re==None or self.Re==0:
                #self.Top_Itr=sp[7] ; self.Bot_Itr=sp[8]
                self.CD=self.CDp
            try:
                self.CL_CD=float(self.CL)/float(self.CD)
            except:
                self.isconv=False
        if(self.isconv==False):
            print("WARNING: unconverged Xfoil run")

    def get_distrib(self):
        ''' Gets the coordinates that were used by xfoil, and the corresponding
        values of pressure coefficient. Stores the results in xout,yout,Cp'''
        self.xout=np.loadtxt(self.cpfile,usecols=(0,))
        self.Cp=np.loadtxt(self.cpfile,usecols=(1,)) 
        # y not stored in previous file, so we retrieve it from extended info
        # IMPORTANT: extinfo file was found out to be "fixed column width" (1+9+9+9)
        yout=[]
        with open(self.extinfo, 'r') as file_:
            for line in file_:
                if(line[0]=="#"):
                    continue
                ystring = line[19:27+1]
                yout.append(float(ystring))
        self.yout=np.array(yout)
    def clean_output(self):
        ''' clean output files and BL files -if there- '''
        for n in [self.coeffs,self.cpfile,self.extinfo,":00.bl","500.bl"]:
            if(os.path.isfile(n)):
                os.remove(n)
    def _write_geometry_file(self):
        ''' Writes airfoil coordinates in temp file'''
        np.savetxt(self.geom, np.c_[self.xin,self.yin])
    def _write_input_file(self):
        ''' Writes the xfoil commands to be executed'''
        f=open(self.rundef,"w")
        # Deactivate graphical interface (necessary when running Ubuntu for Windows or Manjaro, but not native Ubuntu)
        if(self.nograph==True):
            f.write('plop\n')
            f.write('G\n')
            f.write('\n')
        f.write('load %s\n' %(self.geom))
        f.write('\n')
        if self.nrepanel !=None:
            f.write('ppar\n')
            f.write('N\n')
            f.write('%s\n' %(self.nrepanel))
            f.write('\n')
            f.write('\n')   
        f.write('oper\n')
        # Toggle compressibility
        if(self.Mach!=None):
            f.write('Mach %s\n' %(self.Mach))
        # Toggle viscosity
        if(self.Re!=None):
            f.write('visc\n')
            f.write('%s\n' %(self.Re))
        # Transition params
        if(self.n_crit!=None):
            f.write('vpar\n')
            f.write('N %s\n' %(self.n_crit))
            f.write('\n')
        if(self.xtr!=None):
            f.write('vpar\n')
            f.write('xtr\n')
            f.write('%s\n' %(self.xtr))
            f.write('%s\n' %(self.xtr))
            f.write('\n')
        f.write('pacc\n')
        f.write('\n')
        f.write('\n')
        if(self.max_iter!=None):
            f.write('iter %s\n' %(self.max_iter))
        f.write('alfa %s\n' %(self.alphadeg))
        f.write('dump %s\n' %(self.extinfo))
        f.write('cpwr %s\n' %(self.cpfile))
        f.write('pwrt\n')
        f.write('%s\n' %(self.coeffs))
        f.write('\n')
        f.write('quit\n')
        f.close()
    def _run_code(self):
        ''' Executes xfoil '''
        mycmd="xfoil < %s" %(self.rundef)
        os.system(mycmd)
    def _clean_input(self):
        ''' clean input files '''
        for n in [self.rundef]:
            if(os.path.isfile(n)):
                os.remove(n)
