import os
import KratosMultiphysics.kratos_utilities

def Clear():
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('ROM_Snapshots')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('ROM_Skin_Data')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('HROM_Snapshots')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('HHROM_Snapshots')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('HROM_Skin_Data')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('HHROM_Skin_Data')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('RBF_Snapshots')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('RBF_Skin_Data')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Results')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Train_Captures')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Test_Captures')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('case_data.xlsx')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('upwind_elements_list.txt')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('trailing_edge_element_id.txt')
    os.mkdir('ROM_Snapshots')
    os.mkdir('ROM_Skin_Data')
    os.mkdir('HROM_Snapshots')
    os.mkdir('HHROM_Snapshots')
    os.mkdir('HROM_Skin_Data')
    os.mkdir('HHROM_Skin_Data')
    os.mkdir('RBF_Snapshots')
    os.mkdir('RBF_Skin_Data')
    os.mkdir('Train_Captures')
    os.mkdir('Test_Captures')

def ClearAll():
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('FOM_Snapshots')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('FOM_Skin_Data')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('ROM_Snapshots')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('ROM_Skin_Data')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('HROM_Snapshots')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('HHROM_Snapshots')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('HROM_Skin_Data')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('HHROM_Skin_Data')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('RBF_Snapshots')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('RBF_Skin_Data')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Results')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Train_Captures')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Test_Captures')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('rom_data')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('__pycache__')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Clustering_mu_test_plots')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Mu_by_clusters')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('rom_bases')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('FOM_data.xlsx')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('case_data.xlsx')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_train.npy')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_train_list.npy')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_train_not_scaled.npy')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_test.npy')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_test_list.npy')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_test_not_scaled.npy')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('MuValues.png')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('MuValuesNotScaled.png')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('Rom_transonic_potential_flow_Naca0012_clusters.post.lst')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('upwind_elements_list.txt')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('trailing_edge_element_id.txt')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('FOM_SnapshotsMatrix_conv.npy')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('PEBL clustering.png')

if __name__ == "__main__":

    ClearAll()