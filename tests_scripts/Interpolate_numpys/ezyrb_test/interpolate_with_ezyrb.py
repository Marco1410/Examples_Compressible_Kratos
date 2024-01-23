import os
import numpy as np
np.float=float
from ezyrb import Database
from ezyrb import ReducedOrderModel as ROM
from ezyrb import RBF, POD

#----------- (Reading all files with specific extension) ------------
Train_Directory_Path = './solutions_train/snapshots'
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
            parameters_list.append([filename_without_extension.split("_")[-1]])

snapshots_list = np.load("solutions_train/fom_snapshots.npy")

parameters = np.array(parameters_list).reshape(-1,1)
snapshots = np.array(snapshots_list).transpose()

#print("Files' names --> ", filename_list)
print("Number of snapshots = ", count)
print("Parameters shape = ", parameters.shape)
print("Snapshots shape = ", snapshots.shape)
#====================================================================

#-------------------- (Solve for new parameters) --------------------
new_parameters_list = parameters_list #[[0.0,0.7],[0.0,0.8],[0.918,0.7],[2.0,0.7],[2.0,0.8]];
solutions_list = [[0.798],[0.721],[0.754],[0.732],[0.787]] #[[0.338,0.798],[0.838,0.721],[1.338,0.754],[1.588,0.732],[1.838,0.787]];
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
