import os
import KratosMultiphysics.kratos_utilities

def Clear():
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('ROM_Snapshots')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('HROM_Snapshots')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('RBF_Snapshots')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Results')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('case_data.xlsx')
    os.mkdir('ROM_Snapshots')
    os.mkdir('HROM_Snapshots')
    os.mkdir('RBF_Snapshots')

def ClearAll():
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('FOM_Snapshots')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('FOM_Snapshots_nc')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('ROM_Snapshots')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('HROM_Snapshots')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('RBF_Snapshots')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Results')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('rom_data')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('__pycache__')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('FOM_data.xlsx')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('case_data.xlsx')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_train.npy')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_train_list.npy')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_train_not_scaled.npy')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_train_not_scaled_aux.npy')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_test.npy')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_test_list.npy')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_test_not_scaled.npy')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_test_not_scaled_aux.npy')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('MuValues.png')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('MuValuesNotScaled.png')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('ROM_01.post.lst')

if __name__ == "__main__":

    ClearAll()