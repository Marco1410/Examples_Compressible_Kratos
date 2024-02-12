import os,shutil,stat,subprocess

######################################################
# CLASS FOR HANDLING COMMUNICATION IN LOCAL RUNS #####
######################################################

class local:
    """Manages the communication when running simulations locally"""

    def run_command(self,mycommand,runfold=None,envdict=None):
        '''Runs the given command with Linux OS. Returns terminal info.
        If runfold!=None, executes the command in that folder.
        If envdict=None, includes additional content to environment variables
           that already existed.
        NOTE: will not leave the function till the process is done.'''
        # If additional environment variables are needed, we include them
        my_env = os.environ.copy()
        if(envdict!=None):
            for keynow in list(envdict.keys()):
                if keynow in my_env: 
                    my_env[keynow] = "%s:%s" %(my_env[keynow],envdict[keynow])
                else:
                    my_env[keynow] = "%s" %(envdict[keynow])                
        if(runfold==None):
            result = subprocess.run(mycommand, stdout=subprocess.PIPE, 
                env=my_env)
        else:
            result = subprocess.run(mycommand, stdout=subprocess.PIPE,
                env=my_env,cwd=runfold)
        return result.stdout.decode('utf-8')

    def put_data(self,src,dst,isoverwrite=True):
        '''Copies a file/folder from one location (src) to the other (dst).,
           If isoverwrite==False, existing files will be kept'''
        iscopy=True
        if os.path.exists(dst):
            if(isoverwrite==False):
                print("WARNING: %s already exists, keeping old files" %(dst))
                return
        isfold=os.path.isdir(src)
        if isfold:
            shutil.copytree(src,dst)
        else:
            shutil.copy(src,dst)

    def remove_local_data(self,localsrc):
        '''Removes a file/folder from local file system'''
        isfold=os.path.isdir(localsrc)
        if isfold:
            shutil.rmtree(localsrc)
        else:
            os.remove(localsrc)

    def make_localdir_p(self,locfold):
        '''Creates a local directory if it does not exist. Otherwise it does nothing'''
        dirs_ = []
        while len(locfold) > 1:
            dirs_.append(locfold)
            locfold, _  = os.path.split(locfold) 
        if len(locfold) == 1 and not locfold.startswith("/"): 
            dirs_.append(locfold) 
        while len(dirs_):
            dir_ = dirs_.pop()
            if not os.path.isdir(dir_):
                print ("INFO: Making directory %s" %(dir_))
                os.mkdir(dir_)