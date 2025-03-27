import KratosMultiphysics.kratos_utilities

def ClearAll():
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Mu_parameters')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('rom_data')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('__pycache__')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Results')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('output_plots')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Selected_elements_lists')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('FOM_data.xlsx')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('ROM_data.xlsx')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('resume.xlsx')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('selected_elements_list.txt')
    KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('test.post.lst')

if __name__ == "__main__":

    ClearAll()