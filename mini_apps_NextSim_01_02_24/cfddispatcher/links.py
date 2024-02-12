from collections import OrderedDict
from string import Template
from shutil import copytree ,copy
from copy import deepcopy
import os,time
from . import comm

################################################################
## CLASS FOR GROUPS OF SIMULATIONS #############################
################################################################

class simgroup:
    """Common class for groups of simulations"""    
    def add_simulation(self,sim,tag):
        '''Includes a new simulation in the simgroup'''
        if(tag in list(self.simlst.keys())):
            print("WARNING: %s already existed in the simgroup. Overwritten" %(tag))
        self.simlst[tag]=sim

    def remove_simulation(self,tag):
        '''Removes a simulation named as tag from the simlst'''
        if(tag not in list(self.simlst.keys())):
            print("WARNING: %s was not in the simgroup. Not deleted" %(tag))
        self.simlst.pop(tag)

    def deploy_sims_locally(self):
        '''Creates the simulations structure in the local file system,
        and includes the files of the templated folder and the launching file'''
        for simtag in list(self.simlst.keys()):
            print("INFO: Creating simulation %s localy" %(simtag))
            dstfold=os.path.join(self.locfold,simtag)
            # Create simulation folders in the local file system 
            self.comm.make_localdir_p(dstfold)
            thissym=self.simlst[simtag]
            # We first copy the whole tree, and then apply templating one by one
            # NOTE: otherwise, it is hard to deal with subfolders 
            for firstlev in os.listdir(thissym.tfold):
                thisentry=os.path.join(thissym.tfold,firstlev)
                if os.path.isdir(thisentry):
                    # NOTE: We set dirs_exist_ok to overwrite the files, otherwise
                    #       it will simply raise an error, instead of keeping the origin
                    copytree(thisentry,os.path.join(dstfold,firstlev),dirs_exist_ok=True)
                else:
                    copy(thisentry,dstfold)
            # Copying all files included in the template folder
            for firstlev in os.listdir(thissym.tfold):
                subfold=os.path.join(thissym.tfold,firstlev)
                if os.path.isdir(subfold):
                    for secondlev in os.listdir(subfold):
                        if os.path.isdir(secondlev):
                            print("ERROR: The template folder cannot contain 1 level of subfolders")
                            exit()
                        else:
                            srcpath=os.path.join(subfold,secondlev)
                            dstpath=os.path.join(dstfold,firstlev)
                            dstpath=os.path.join(dstpath,secondlev)
                            apply_template(srcpath,dstpath,thissym.tdict)
                else:
                    dstpath=os.path.join(dstfold,firstlev)
                    apply_template(subfold,dstpath,thissym.tdict)
            # Create launching file
            thissym.launcher.write_launching_file(dstfold)

###############################################################
# CLASS FOR HANDLING A GROUP OF LOCAL SIMULATIONS ############
###############################################################

class localsimgroup(simgroup):
    """Handles group of local simulations"""

    def __init__(self,locfold):
        '''Class initialization. Builds a local comm class,
        and relates the simulation to a local folder (locfold)'''
        self.comm=comm.local()
        self.locfold=locfold
        # We create the destination folder if it does not exist
        self.comm.make_localdir_p(locfold)
        # Lists of simulations included in this group
        self.simlst=OrderedDict()

    def perform_simulations(self,envdict=None):
        '''Tells the cluster to run every single simulation of the group.
        Use envdict!=None to include more content in already existing environment variables'''
        # We prepare a dict to track the status of the runs
        self.jobids=OrderedDict()
        for simtag in list(self.simlst.keys()):
            dstfold=os.path.join(self.locfold,simtag)
            thissym=self.simlst[simtag]
            print("INFO: Performing simulation %s locally" %(simtag))
            jobid=thissym.launcher.perform_computation(comm=self.comm,
                executefold=dstfold,envdict=envdict)

    def copy_additional_data(self,keepsourcedata=True):
        '''Copies all the additional data that is necessary for the simulation, besides
        the templated files. The mesh is typically one of such files.
        Set keepsourcedata=False if you want to remove the files from their initial
        locations'''
        for simtag in list(self.simlst.keys()):
            print("INFO: Copying additional files of %s" %(simtag))
            dstfold=os.path.join(self.locfold,simtag)
            thissym=self.simlst[simtag]
            # Loop in additional input data
            for srcpath in thissym.inputdata:
                dataname=os.path.basename(srcpath)
                dstpath=os.path.join(dstfold,dataname)
                self.comm.put_data(srcpath,dstpath)
        # We delete the data after all the copies has been done,
        # in case there were files shared by several computations
        if not keepsourcedata:
            print("INFO: Deleting additional simulation files from initial locations")
            for simtag in list(self.simlst.keys()):
                thissym=self.simlst[simtag]
                # Loop in additional data
                for datapath in thissym.inputdata:
                    if(os.path.exists(datapath)):
                        self.comm.remove_local_data(datapath)

##############################################
### A SINGLE SIMULATION ######################
############################################## 

class sim:
    """Class of a single computation. tfold is the template folder for this simulation,
    and tdict contains the variables to replace"""
    def __init__(self,tfold,tdict):
        '''Class initialization'''
        self.tfold=tfold
        # Complement tdict with main cfdispatcher folder to allow importing it
        # IMPORTANT: Then the variable $cfddfold is available in all templates to be imported
        tdictext=deepcopy(tdict)
        thisfold=os.path.dirname(os.path.realpath(__file__))
        parentfold=os.path.abspath(os.path.join(thisfold, os.pardir))
        tdictext["cfddfold"]=parentfold
        self.tdict=tdictext
        # Input files and folders         
        self.inputdata=[]

    def add_input_data(self,datan):
        '''Adds a file/folder to the input of the simulation'''
        self.inputdata.append(datan)

    def add_launcher(self,launcher):
        '''Relates a simulation to a launcher (from queue.py),
        that could be for instance a SLURM class'''
        self.launcher=launcher

##############################################
### AUXILIAR FUNCTIONS 
##############################################

def apply_template(srcpath,dstpath,dictsubs):
    '''applies the template of dictsubs to srcpath, and stores it in dstpath'''
    fin=open(srcpath,'r')
    data=fin.read()
    fin.close()
    t = Template(data)
    s=t.substitute(dictsubs)
    # Copy replaced content
    fout=open(dstpath,'w')
    fout.write(s)
    fout.close() 
