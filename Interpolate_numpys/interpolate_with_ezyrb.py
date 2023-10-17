import os
import numpy as np
np.float=float
from ezyrb import Database
from ezyrb import ReducedOrderModel as ROM
from ezyrb import RBF, POD

#----------- (Reading all files with specific extension) ------------
Train_Directory_Path = './solutions_train'
Test_Directory_Path = './solutions_test'
Files_Extension = '.npy'

count = 0
filename_list = []; parameters_list = []; snapshots_list = []
for path, folders, files in os.walk(Train_Directory_Path):
    for filename in sorted(files): # Sorting files by name
        filename_without_extension, extension = os.path.splitext(filename)
        if extension == Files_Extension:
            count +=1
            full_path = '{}/{}'.format(path,filename)
            filename_list.append(filename)
            parameters_list.append(int(filename_without_extension.split("_")[-1]))
            snapshots_list.append(np.load(full_path))

parameters = np.array(parameters_list).reshape(-1,1)
snapshots = np.array(snapshots_list)

print("Files' names --> ", filename_list)
print("Number of snapshots = ", count)
print("Parameters shape = ", parameters.shape)
print("Snapshots shape = ", snapshots.shape)
#====================================================================

#-------------------- (Solve for new parameters) --------------------
new_parameters_list = [45,55,65,75,85,95];
solutions_list = [10,20,30,40,50,60];
db = Database(parameters, snapshots)
pod = POD()
rbf = RBF()
rom = ROM(db, pod, rbf)
rom.fit()
for element in solutions_list:
    norm_error_absolut = np.linalg.norm(np.load(Train_Directory_Path+f'/trainStepSolution_{element}.npy')-rom.predict([element]))
    norm_original_solution =  np.linalg.norm(np.load(Train_Directory_Path+f'/trainStepSolution_{element}.npy'))
    print('train error: ', norm_error_absolut/norm_original_solution*100, ' %' )
for element in new_parameters_list:
    norm_error_absolut = np.linalg.norm(np.load(Test_Directory_Path+f'/testStepSolution_{element}.npy')-rom.predict([element]))
    norm_original_solution =  np.linalg.norm(np.load(Test_Directory_Path+f'/testStepSolution_{element}.npy'))
    print('test error: ', norm_error_absolut/norm_original_solution*100, ' %' )
#====================================================================
