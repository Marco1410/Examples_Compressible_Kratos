import KratosMultiphysics
import numpy as np
from scipy.stats import qmc
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import pickle
import os
import KratosMultiphysics.kratos_utilities
from KratosMultiphysics.MeshMovingApplication.mesh_moving_analysis import MeshMovingAnalysis
from KratosMultiphysics.RomApplication.rom_manager import RomManager

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def GenerateMeshes(shape_parameters,test_shapes,plot_airfoils):
    if plot_airfoils:
        fig = plt.figure()
        fig.subplots_adjust(top=0.8)
        fig.set_figwidth(10.0)
        fig.set_figheight(7.0)

    for i in range(len(shape_parameters)):    
        with open("ProjectParametersMeshAirfoil.json",'r') as parameter_file:
            airfoil_mesh_parameters = KratosMultiphysics.Parameters(parameter_file.read()) 
        airfoil_mesh_parameters["output_processes"].RemoveValue("gid_output")
        airfoil_mesh_parameters["solver_settings"]["model_import_settings"]["input_type"].SetString("mdpa")
        airfoil_mesh_parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString("naca0012_0aoa_00")
        
        airfoil_model = KratosMultiphysics.Model()
        airfoil_mesh_simulation = MeshMovingAnalysis(airfoil_model,airfoil_mesh_parameters)
        airfoil_mesh_simulation.Initialize()

        modelpart = airfoil_model["MainModelPart.Body2D_Body"]
        if plot_airfoils:
            xp = np.zeros(modelpart.NumberOfNodes())
            yp = np.zeros(modelpart.NumberOfNodes())
        for id,node in enumerate(modelpart.Nodes):
            X,Y,Z = coordinates_new_shape(node.X0,node.Y0,node.Z0,shape_parameters[i,:])
            node.SetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT_X,X-node.X0)
            node.SetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT_Y,Y-node.Y0)
            node.SetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT_Z,Z-node.Z0)
            node.Fix(KratosMultiphysics.MESH_DISPLACEMENT_X)
            node.Fix(KratosMultiphysics.MESH_DISPLACEMENT_Y)
            node.Fix(KratosMultiphysics.MESH_DISPLACEMENT_Z)

            if plot_airfoils:
                xp[id-1] = X
                yp[id-1] = Y

        airfoil_mesh_simulation.RunSolutionLoop()
        airfoil_mesh_simulation.Finalize()

        mainmodelpart = airfoil_model["MainModelPart"]
        if (i>(len(shape_parameters)-test_shapes-1)):
            id_test = test_shapes - (len(shape_parameters)-i)
            KratosMultiphysics.ModelPartIO("Meshes/test_mesh_"+str(id_test), KratosMultiphysics.IO.WRITE | KratosMultiphysics.IO.MESH_ONLY).WriteModelPart(mainmodelpart)
            fig = plt.plot( xp, yp, "+", markersize=5, label = "Test  Airfoil "+str(id_test)+" t="+str(np.round(shape_parameters[i,0],3))+" m="+str(np.round(shape_parameters[i,1],3))+" p="+str(np.round(shape_parameters[i,2],3)))
            # guardar aqui las coordenadas de los perfiles
            fout=open("Data/mesh_"+str(id_test)+".dat", 'w')
            for j in range(modelpart.NumberOfNodes()):
                fout.write("%s %s\n" %(xp[j],yp[j]))
            fout.close()
        else:
            fig = plt.plot( xp, yp, "+", markersize=5, label = "Train Airfoil "+str(i)+" t="+str(np.round(shape_parameters[i,0],3))+" m="+str(np.round(shape_parameters[i,1],3))+" p="+str(np.round(shape_parameters[i,2],3)))
            KratosMultiphysics.ModelPartIO("Meshes/mesh_"+str(i), KratosMultiphysics.IO.WRITE | KratosMultiphysics.IO.MESH_ONLY).WriteModelPart(mainmodelpart)
            # guardar aqui las coordenadas de los perfiles            
            fout=open("Data/mesh_"+str(i)+".dat", 'w')
            for j in range(modelpart.NumberOfNodes()):
                fout.write("%s %s\n" %(xp[j],yp[j]))
            fout.close()

    if plot_airfoils:
        fig = plt.axis([-0.1,1.2,-0.3,0.3])
        fig = plt.grid()
        fig = plt.legend()
        fig = plt.title('Airfoils')
        fig = plt.xlabel("x")
        fig = plt.ylabel("y")
        fig = plt.tight_layout()
        fig = plt.savefig("Images/airfoils.png")

def coordinates_new_shape(x,y,z,shape_parameters):
    # Definir el perfil base: Naca Serie 4 digitos
    # Longitud de la cuerda c = 1
    # n Número de puntos en el perfil
    # x vector de posiciones a lo largo de la cuerda
    # Definir el rango de valores para el espesor máximo y curvatura
    # naca xxxx posición del espesor máximo = 30%
    # Ejemplo
    # tmax_range = [0.10, 0.15] # Espesor máximo 
    # m_range    = [0.0 , 0.03] # Máxima curvatura
    # p_range    = [0.25, 0.35] # Posicion de máxima curvatura
    tmax,m,p = shape_parameters 
    if x < 0.0:
        x = -x
    # Calcular linea de curvatura actual
    if x <= p:
        yc = (m/p**2)*(2*p*x-x**2)
        dyc = (2*m/p**2)*(p-x)
    else: 
        yc = (m/(1-p)**2)*((1-2*p)+2*p*x-x**2)
        dyc = (2*m/(1-p)**2)*(p-x)
    # Calcular el espesor máximo para el perfil actual
    yt = 5 * tmax * (0.2969*np.sqrt(x) - 0.1260*(x) - 0.3516*(x)**2 + 0.2843*(x)**3 - 0.1036*(x)**4)
    # Calcular las coordenadas del perfil actual
    if y > 0:
        xp = x - yt * np.sin(np.arctan(dyc)) 
        yp = yc + yt * np.cos(np.arctan(dyc))
    else:
        xp = x + yt * np.sin(np.arctan(dyc))
        yp = yc - yt * np.cos(np.arctan(dyc))
    return xp, yp, z

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def CustomizeSimulation(cls, global_model, parameters):

    class CustomSimulation(cls):

        def __init__(self, model,project_parameters, custom_param = None):
            super().__init__(model,project_parameters)
            self.custom_param  = custom_param

        def Initialize(self):

            angle_of_attack = parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].GetDouble()
            parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(0.0) 
            with open("ProjectParametersMeshMoving.json",'r') as parameter_file:
                mesh_parameters = KratosMultiphysics.Parameters(parameter_file.read()) 
            mesh_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["rotation_angle"].SetDouble(angle_of_attack)
            mesh_parameters["output_processes"].RemoveValue("gid_output")
            mesh_simulation = MeshMovingAnalysis(self.model,mesh_parameters)

            super().Initialize()
            
            mesh_simulation.Run()

        def FinalizeSolutionStep(self):
            super().FinalizeSolutionStep()
            # guardar aqui datos directamente de la skin y plotea
            id_case = np.loadtxt("Data/case.dat",usecols=(0,))
            mesh_id = np.loadtxt("Data/case.dat",usecols=(1,))
            test = np.loadtxt("Data/case.dat",usecols=(2,))
            name = np.loadtxt("Data/case.dat",usecols=(3,))
            prefix = np.loadtxt("Data/case.dat",usecols=(4,))
            if prefix == 0:
                if name == 1:
                    case = "FOM"
                if name == 2:
                    case = "ROM"
                if name == 3:
                    case = "HROM"
                KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('Data/case')
                if test == 1:
                    fout=open("Data/test_mesh_"+str(int(mesh_id))+"_"+case+"_"+str(int(id_case))+".dat",'w')
                else:
                    fout=open("Data/mesh_"+str(int(mesh_id))+"_"+case+"_"+str(int(id_case))+".dat",'w')
                modelpart = self.model["MainModelPart.Body2D_Body"]
                for node in modelpart.Nodes:
                    x=node.X ; y=node.Y ; z=node.Z
                    cp=node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
                    fout.write("%s %s %s %s\n" %(x,y,z,cp))
                fout.close()
            else:
                if name == 1:
                    case = "FOM"
                if name == 2:
                    case = "ROM"
                if name == 3:
                    case = "HROM"
                KratosMultiphysics.kratos_utilities.DeleteFileIfExisting('Data/case')
                if test == 1:
                    fout=open("Data/test_mesh_"+str(int(mesh_id))+"_"+"Adjoint"+"_"+case+"_"+str(int(id_case))+".dat",'w')
                else:
                    fout=open("Data/mesh_"+str(int(mesh_id))+"_"+"Adjoint"+"_"+case+"_"+str(int(id_case))+".dat",'w')
                modelpart = self.model["MainModelPart.Body2D_Body"]
                for node in modelpart.Nodes:
                    x=node.X ; y=node.Y ; z=node.Z
                    ssx=node.GetSolutionStepValue(KratosMultiphysics.SHAPE_SENSITIVITY_X)
                    ssy=node.GetSolutionStepValue(KratosMultiphysics.SHAPE_SENSITIVITY_Y)
                    ssz=node.GetSolutionStepValue(KratosMultiphysics.SHAPE_SENSITIVITY_Z)
                    fout.write("%s %s %s %s %s %s\n" %(x,y,z,ssx,ssy,ssz))
                fout.close()

        def CustomMethod(self):
            return self.custom_param

    return CustomSimulation(global_model, parameters)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def UpdateProjectParameters(parameters, mu=None, Id=None, name=None, mesh_id=None, test=None, prefix=None):
    if prefix == "Primal_":
        prefix = 0
    else:
        prefix = 1        
    if name == "FOM" or name == "FOM_test":
        case = 1
    if name == "ROM" or name == "ROM_test":
        case = 2
    if name == "HROM" or name == "HROM_test":
        case = 3
    if test:
        parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString("Meshes/test_mesh_"+str(mesh_id))
        fout=open("Data/case.dat",'w')
        fout.write("%s %s %s %s %s\n" %(Id,mesh_id,1,case,prefix))
        fout.close()
    else:
        parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString("Meshes/mesh_"+str(mesh_id))
        fout=open("Data/case.dat",'w')
        fout.write("%s %s %s %s %s\n" %(Id,mesh_id,0,case,prefix))
        fout.close()
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(mu[0])
    parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mu[1])
    parameters["processes"]["list_other_processes"][0]["Parameters"]["file_settings"]["file_name"].SetString(name+"_Primal_Data/MainModelPart_"+str(Id))
    return parameters

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def GetRomManagerParameters():
    general_rom_manager_parameters = KratosMultiphysics.Parameters("""{
            "rom_stages_to_train" : ["ROM","HROM"],      // ["ROM","HROM"]
            "rom_stages_to_test"  : ["ROM","HROM"],      // ["ROM","HROM"]
            "paralellism" : null,                        // null, TODO: add "compss"
            "projection_strategy": "galerkin",           // "lspg", "galerkin", "petrov_galerkin"
            "save_gid_output": false,                     // false, true #if true, it must exits previously in the ProjectParameters.json
            "save_vtk_output": false,                    // false, true #if true, it must exits previously in the ProjectParameters.json
            "output_name": "id",                         // "id" , "mu"
            "ROM":{
                "svd_truncation_tolerance": 1e-9,
                "model_part_name": "MainModelPart",                                      // This changes depending on the simulation: Structure, FluidModelPart, ThermalPart #TODO: Idenfity it automatically
                "nodal_unknowns": [ "VELOCITY_POTENTIAL",
                                    "AUXILIARY_VELOCITY_POTENTIAL"],                     // Main unknowns. Snapshots are taken from these
                "rom_basis_output_format": "json",                                       // "json" "numpy"
                "rom_basis_output_name": "ROM_Basis/PrimalRomParameters",
                "snapshots_control_type": "step",                                        // "step", "time"
                "snapshots_interval": 1,
                "petrov_galerkin_training_parameters":{
                    "basis_strategy": "residuals",                                        // 'residuals', 'jacobian'
                    "include_phi": false,
                    "svd_truncation_tolerance": 1e-9,
                    "echo_level": 0
                }
            },
            "HROM":{
                "element_selection_type": "empirical_cubature",
                "element_selection_svd_truncation_tolerance": 1e-12,
                "create_hrom_visualization_model_part" : false,
                "echo_level" : 0
            }
        }""")
    return general_rom_manager_parameters

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# mult params
#     
def get_multiple_params_by_Halton_test(number_of_values,angle,mach):
    sampler = qmc.Halton(d=2)
    sample = sampler.random(number_of_values)
    #Angle of attack signo segun el eje de rotacion
    mu = []
    values = qmc.scale(sample, [angle[0],mach[0]], [angle[1],mach[1]])
    for i in range(number_of_values):
        #Angle of attack , Mach infinit
        mu.append([values[i,0] * np.pi / 180.0, values[i,1]])
    return mu

def get_multiple_params_by_Halton_train(number_of_values,angle,mach):
    sampler = qmc.Halton(d=2)
    sample = sampler.random(number_of_values)
    mu = []
    values = qmc.scale(sample, [angle[0],mach[0]], [angle[1],mach[1]])
    values[0,0] = angle[0]
    values[0,1] = mach[0]
    values[1,0] = angle[0]
    values[1,1] = mach[1]
    values[number_of_values-1,0] = angle[1]
    values[number_of_values-1,1] = mach[1]
    values[number_of_values-2,0] = angle[1]
    values[number_of_values-2,1] = mach[0]
    for i in range(number_of_values):
        #Angle of attack , Mach infinit
        mu.append([values[i,0] * np.pi / 180.0, values[i,1]])
    return mu

def get_multiple_shape_params_by_Halton(number_of_values,tmax,m,p):
    sampler = qmc.Halton(d=3)
    sample = sampler.random(number_of_values)
    return qmc.scale(sample,[tmax[0],m[0],p[0]],[tmax[1],m[1],p[1]])
    
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# plot
# 
def plot_mu_values(mu_train,mu_test,xl,xu,yl,yu):
    mu_train_a = np.zeros(len(mu_train))
    mu_train_m = np.zeros(len(mu_train))
    mu_test_a  = np.zeros(len(mu_test))
    mu_test_m  = np.zeros(len(mu_test))
    fig = plt.figure()
    fig.subplots_adjust(top=0.8)
    fig.set_figwidth(10.0)
    fig.set_figheight(7.0)
    for i in range(len(mu_train)):
        mu_train_a[i] = -mu_train[i][0] * 180 / np.pi
        mu_train_m[i] = mu_train[i][1]
        fig = plt.plot(mu_train_m[i], mu_train_a[i], 'bs', markersize = 5.0, label="Train "+str(np.round(mu_train_a[i],1))+"º "+str(np.round(mu_train_m[i],2)))
    for i in range(len(mu_test)):
        mu_test_a[i] = -mu_test[i][0] * 180 / np.pi
        mu_test_m[i] = mu_test[i][1]
        fig = plt.plot(mu_test_m[i], mu_test_a[i], 'ro', markersize = 5.0, label="Test  "+str(np.round(mu_test_a[i],1))+"º "+str(np.round(mu_test_m[i],2)))
    fig = plt.title('Mu Values')
    fig = plt.axis([xl-0.1,xu+0.2,-yu-0.1,-yl+0.1])
    fig = plt.ylabel('Alpha')
    fig = plt.xlabel('Mach')
    fig = plt.grid()
    fig = plt.legend()
    fig = plt.tight_layout()
    fig = plt.savefig("Images/MuValues.png")

def plot_Cps(mu_train,mu_test,NumberOfTrainShapes,NumberOfTestShapes,NumberofMuTrain,NumberOfMuTest):
    case_names = ["FOM","ROM","HROM"]
    for i in range(NumberOfTrainShapes):
        for j in range(NumberofMuTrain):
            fig = plt.close('all')
            fig = plt.figure()
            fig.subplots_adjust(top=0.8)
            fig.set_figwidth(10.0)
            fig.set_figheight(7.0)
            for name in case_names:
                x  = np.loadtxt("Data/mesh_"+str(i)+"_"+name+"_"+str(j)+".dat",usecols=(0,))
                cp = np.loadtxt("Data/mesh_"+str(i)+"_"+name+"_"+str(j)+".dat",usecols=(3,))
                fig = plt.plot(x, cp, '+', markersize = 3.0, label=name+" "+str(np.round(-mu_train[j][0] * 180 / np.pi,1))+"º " + str(np.round(mu_train[j][1],2)))
            fig = plt.title('Cp vs x')
            fig = plt.axis([-0.1,1.2,1.1,-3.0])
            fig = plt.ylabel('Cp')
            fig = plt.xlabel('x')
            fig = plt.grid()
            fig = plt.legend()
            fig = plt.tight_layout()
            fig = plt.savefig("Images/mesh_"+str(i)+"_"+str(j)+".png")

    for i in range(NumberOfTestShapes):
        for j in range(NumberOfMuTest):
            fig = plt.close('all')
            fig = plt.figure()
            fig.subplots_adjust(top=0.8)
            fig.set_figwidth(10.0)
            fig.set_figheight(7.0)
            for name in case_names:
                x  = np.loadtxt("Data/test_mesh_"+str(i)+"_"+name+"_"+str(j)+".dat",usecols=(0,))
                cp = np.loadtxt("Data/test_mesh_"+str(i)+"_"+name+"_"+str(j)+".dat",usecols=(3,))
                fig = plt.plot(x, cp, '+', markersize = 3.0, label=name+" "+str(np.round(-mu_test[j][0] * 180 / np.pi,1))+"º " + str(np.round(mu_test[j][1],2)))
            fig = plt.title('Cp vs x')
            fig = plt.axis([-0.1,1.2,1.1,-3.0])
            fig = plt.ylabel('Cp')
            fig = plt.xlabel('x')
            fig = plt.grid()
            fig = plt.legend()
            fig = plt.tight_layout()
            fig = plt.savefig("Images/test_mesh_"+str(i)+"_"+str(j)+".png")

def plot_Sensitivities(mu_train,mu_test,NumberOfTrainShapes,NumberOfTestShapes,NumberofMuTrain,NumberOfMuTest):
    case_names = ["FOM","ROM","HROM"]
    for i in range(NumberOfTrainShapes):
        for j in range(NumberofMuTrain):
            fig = plt.close('all')
            fig = plt.figure()
            fig.subplots_adjust(top=0.8)
            fig.set_figwidth(10.0)
            fig.set_figheight(7.0)
            for name in case_names:
                x  = np.loadtxt("Data/mesh_"+str(i)+"_"+"Adjoint"+"_"+name+"_"+str(j)+".dat",usecols=(0,))
                ssx = np.loadtxt("Data/mesh_"+str(i)+"_"+"Adjoint"+"_"+name+"_"+str(j)+".dat",usecols=(3,))
                ssy = np.loadtxt("Data/mesh_"+str(i)+"_"+"Adjoint"+"_"+name+"_"+str(j)+".dat",usecols=(4,))
                ssz = np.loadtxt("Data/mesh_"+str(i)+"_"+"Adjoint"+"_"+name+"_"+str(j)+".dat",usecols=(5,))
                fig = plt.plot(x, np.sqrt(ssx**2+ssy**2+ssz**2), '+', markersize = 3.0, label=name+" "+str(np.round(-mu_train[j][0] * 180 / np.pi,1))+"º " + str(np.round(mu_train[j][1],2)))
            fig = plt.title('Shape Sentivitity Norm vs x')
            fig = plt.axis([-0.1,1.2,-0.1,2.0])
            fig = plt.ylabel('Shape Sentivitity Norm')
            fig = plt.xlabel('x')
            fig = plt.grid()
            fig = plt.legend()
            fig = plt.tight_layout()
            fig = plt.savefig("Images/mesh_"+"Adjoint"+"_"+str(i)+"_"+str(j)+".png")

    for i in range(NumberOfTestShapes):
        for j in range(NumberOfMuTest):
            fig = plt.close('all')
            fig = plt.figure()
            fig.subplots_adjust(top=0.8)
            fig.set_figwidth(10.0)
            fig.set_figheight(7.0)
            for name in case_names:
                x  = np.loadtxt("Data/test_mesh_"+str(i)+"_"+"Adjoint"+"_"+name+"_"+str(j)+".dat",usecols=(0,))
                ssx = np.loadtxt("Data/mesh_"+str(i)+"_"+"Adjoint"+"_"+name+"_"+str(j)+".dat",usecols=(3,))
                ssy = np.loadtxt("Data/mesh_"+str(i)+"_"+"Adjoint"+"_"+name+"_"+str(j)+".dat",usecols=(4,))
                ssz = np.loadtxt("Data/mesh_"+str(i)+"_"+"Adjoint"+"_"+name+"_"+str(j)+".dat",usecols=(5,))
                fig = plt.plot(x, np.sqrt(ssx**2+ssy**2+ssz**2), '+', markersize = 3.0, label=name+" "+str(np.round(-mu_test[j][0] * 180 / np.pi,1))+"º " + str(np.round(mu_test[j][1],2)))
            fig = plt.title('Shape Sentivitity Norm vs x')
            fig = plt.axis([-0.1,1.2,-0.1,2.0])
            fig = plt.ylabel('Shape Sentivitity Norm')
            fig = plt.xlabel('x')
            fig = plt.grid()
            fig = plt.legend()
            fig = plt.tight_layout()
            fig = plt.savefig("Images/test_mesh_"+"Adjoint"+"_"+str(i)+"_"+str(j)+".png")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# save / load
# 
def save_mu_parameters(mu_train,mu_test):
    archivo = open('Data/mu_train.dat', 'wb')
    pickle.dump(mu_train, archivo)
    archivo.close()
    archivo = open('Data/mu_test.dat', 'wb')
    pickle.dump(mu_test, archivo)
    archivo.close()

def load_mu_parameters():
    archivo = open('Data/mu_train.dat', 'rb')
    lista = pickle.load(archivo)
    mu_train = np.asarray(lista)
    archivo.close()
    archivo = open('Data/mu_test.dat', 'rb')
    lista = pickle.load(archivo)
    mu_test = np.asarray(lista)
    archivo.close()
    return mu_train,mu_test
    
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
def Clean():
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('PrimalResults')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('AdjointResults')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Data')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('FOM_Primal_Data')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('ROM_Primal_Data')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('HROM_Primal_Data')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Images')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('Meshes')
    KratosMultiphysics.kratos_utilities.DeleteDirectoryIfExisting('ROM_Basis')
    os.mkdir("PrimalResults")
    os.mkdir("AdjointResults")
    os.mkdir("Data")
    os.mkdir("FOM_Primal_Data")
    os.mkdir("ROM_Primal_Data")
    os.mkdir("HROM_Primal_Data")
    os.mkdir("Images")
    os.mkdir("Meshes")
    os.mkdir("ROM_Basis")
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

if __name__ == "__main__":

    Clean()

    NumberOfTrainShapes = 2
    NumberOfTestShapes = 1

    # Definir el rango de valores para el espesor máximo y curvatura
    # naca xxxx posición del espesor máximo = 30%
    tmax_range = [0.10, 0.15] # Espesor máximo 
    m_range    = [0.0 , 0.03] # Máxima curvatura
    p_range    = [0.25, 0.35] # Posicion de máxima curvatura
    shape_parameters = get_multiple_shape_params_by_Halton(NumberOfTrainShapes+NumberOfTestShapes,tmax_range,m_range,p_range)
    GenerateMeshes(shape_parameters,NumberOfTestShapes,plot_airfoils=True)

    NumberofMuTrain = 5 # Minimo 5 por las cuatro esquinas y un punto interno
    NumberOfMuTest  = 1

    # Definir rango de valores de mach y angulo de ataque
    mach_range  = [0.03,0.6]
    angle_range = [-6.0,1.0] #signo segun el eje de rotacion
    mu_train = get_multiple_params_by_Halton_train(NumberofMuTrain,angle_range,mach_range) #fija los puntos en las esquinas
    mu_test  = get_multiple_params_by_Halton_test(NumberOfMuTest,angle_range,mach_range) 

    save_mu_parameters(mu_train,mu_test)
    # mu_train, mu_test = load_mu_parameters()

    plot_mu_values(mu_train,mu_test,mach_range[0],mach_range[1],angle_range[0],angle_range[1])


    general_rom_manager_parameters = GetRomManagerParameters()

    primal_project_parameters_name = "ProjectParametersPrimalROM.json"

    primal_rom_manager = RomManager(primal_project_parameters_name,general_rom_manager_parameters,CustomizeSimulation,UpdateProjectParameters,"Primal_")

    primal_rom_manager.Fit(mu_train,NumberOfTrainShapes) 

    primal_rom_manager.Test(mu_test,NumberOfTestShapes) 

    primal_rom_manager.PrintErrors()

    plot_Cps(mu_train,mu_test,NumberOfTrainShapes,NumberOfTestShapes,NumberofMuTrain,NumberOfMuTest)

    input("press any key to continue with Adjoint problem")

    adjoint_project_parameters_name = "ProjectParametersAdjointROM.json"

    general_rom_manager_parameters["ROM"]["nodal_unknowns"].SetStringArray([ "ADJOINT_VELOCITY_POTENTIAL"
                                                                            ,"ADJOINT_AUXILIARY_VELOCITY_POTENTIAL"])

    general_rom_manager_parameters["ROM"]["rom_basis_output_name"].SetString("ROM_Basis/AdjointRomParameters")

    adjoint_rom_manager = RomManager(adjoint_project_parameters_name,general_rom_manager_parameters,CustomizeSimulation,UpdateProjectParameters,"Adjoint_")

    adjoint_rom_manager.Fit(mu_train,NumberOfTrainShapes) 

    adjoint_rom_manager.Test(mu_test,NumberOfTestShapes) 

    adjoint_rom_manager.PrintErrors()

    plot_Sensitivities(mu_train,mu_test,NumberOfTrainShapes,NumberOfTestShapes,NumberofMuTrain,NumberOfMuTest)

