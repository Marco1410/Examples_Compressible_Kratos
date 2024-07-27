import os
import KratosMultiphysics.kratos_utilities

def Clear():
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Results')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Mu_Captures')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Validation')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('case_data.xlsx')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('FOM_cp_snapshots')
    os.mkdir('Mu_Captures')
    os.mkdir('Validation')

def ClearAll():
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('FOM_Snapshots')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('FOM_Skin_Data')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Results')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Mu_Captures')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('rom_data')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Validation')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('__pycache__')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('FOM_cp_snapshots')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('FOM_data.xlsx')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('case_data.xlsx')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu.npy')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('mu_not_scaled.npy')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('MuValues.png')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('MuValuesNotScaled.png')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('naca0012_dabase_generation_mesh_1.post.lst')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('upwind_elements_list.txt')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('trailing_edge_element_id.txt')

if __name__ == "__main__":

    ClearAll()