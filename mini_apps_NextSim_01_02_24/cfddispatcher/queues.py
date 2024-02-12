import os,re,stat

###########################################################
# CLASS TO MANAGE THE SIMULATIONS THAT ARE SENT LOCALLY  ##
###########################################################

class mypc():
    """Class to manage the run of simulations locally"""
    def __init__(self):
        '''Class initialization'''
        self.launchfile="launch.sh"

    def define_run_command(self,instr):
        '''Gets the run command to include in the bash file'''
        self.run_command=instr

    def write_launching_file(self,localfold):
        '''Writes out a launching file in a local folder'''
        data="%s\n%s\n%s" %(self.generate_header(),
                            self.run_command,
                            self.generate_footer())
        srcpath=os.path.join(localfold,self.launchfile)
        if not os.path.exists(srcpath):
            fout=open(srcpath,'w')
            fout.write(data)
            fout.close()

    def generate_header(self):
        '''Generates the header for the launching file'''
        data='#!/bin/bash\n'
        return data

    def generate_footer(self):
        '''Generates the footer for the launching file'''
        return ""

    def perform_computation(self,comm,executefold,envdict=None):
        '''Runs a computation in executefold through a local comm instance.
        Use envdict!=None to add content to existing environment variables
        NOTE: waits for the computation before returning'''
        # Update permissions so that we can execute it
        execpath=os.path.join(executefold,self.launchfile)
        st = os.stat(execpath)
        os.chmod(execpath, st.st_mode | stat.S_IEXEC)
        # Run executable
        torun="./%s" %(self.launchfile)
        msg_=comm.run_command(torun,runfold=executefold,envdict=envdict)