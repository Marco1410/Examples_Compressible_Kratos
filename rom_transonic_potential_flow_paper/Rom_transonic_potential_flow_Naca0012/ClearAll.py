import KratosMultiphysics.kratos_utilities

def ClearAll():
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('FOM_Snapshots')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('FOM_Skin_Data')

    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('ROM_Snapshots')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('ROM_Skin_Data')
    
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('RBF_Snapshots')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('RBF_Skin_Data')
    
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Train_Captures')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Test_Captures')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Validation')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Results')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Modes')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Mu_history')

    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('rom_data')
    
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('FOM_data.xlsx')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('ROM_data.xlsx')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('resume.xlsx')

    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('selected_elements_list.txt')

    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('Mu_parameters_cp.pdf')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('Singular_values_decay.pdf')

    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('Rom_transonic_potential_flow_Naca0012.post.lst')

if __name__ == "__main__":

    ClearAll()