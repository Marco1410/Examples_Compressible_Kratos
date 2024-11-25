import KratosMultiphysics.kratos_utilities

def ClearAll():
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('FOM_Snapshots')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('FOM_Skin_Data')

    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('ROM_Snapshots')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('ROM_Skin_Data')
    
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('HROM_Snapshots')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('HROM_Skin_Data')
    
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('HHROM_Snapshots')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('HHROM_Skin_Data')
    
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('RBF_Snapshots')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('RBF_Skin_Data')
    
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Train_Captures')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Test_Captures')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Validation')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Results')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Modes')

    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('rom_data')

    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('ResidualsSnapshotsMatrix.zarr')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('PhiHROMMatrix.zarr')
    
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('__pycache__')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('dask-scratch-space')

    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('FOM_data.xlsx')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('case_data.xlsx')

    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_train.npy')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_train_new.npy')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_train_not_scaled.npy')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_train_not_scaled_new.npy')

    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_test.npy')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_test_new.npy')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_test_not_scaled.npy')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_test_not_scaled_new.npy')
    
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_validation.npy')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_validation_not_scaled.npy')

    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('MuValues.pdf')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('MuValuesNotScaled.pdf')

    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('upwind_elements_list.txt')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('wake_elements_list.txt')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('kutta_elements_list.txt')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('trailing_edge_elements_list.txt')

    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('Mu_parameters_cp.pdf')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('Singular_values_decay.pdf')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('EmpiricalCubaturePlot.pdf')

    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('Rom_transonic_potential_flow_Naca0012_C1_2.post.lst')

if __name__ == "__main__":

    ClearAll()