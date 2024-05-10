#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.11.0 with dump python functionality
###

import os

path = os.getcwd()
salome_output_path = 'Mesh'
if not os.path.exists(salome_output_path):
    os.makedirs(salome_output_path)

coordinates_file_name = 'airfoilCoordinates/k1.dat'
AOA = 0.0
Airfoil_MeshSize = 1e-5
Domain_Size = 100
Ratio = 1.001
Growth_Rate = 0.05
Domain_Length = Domain_Size
Domain_Width  = Domain_Size

FarField_MeshSize = Domain_Length / 50.0

print ('Domain_Size = ', Domain_Length, 'AOA = ', AOA, 'FarField_MeshSize = ', FarField_MeshSize, 'Airfoil_MeshSize', Airfoil_MeshSize)

import sys
import salome# type: ignore
salome.salome_init()
theStudy = salome.myStudy
import salome_notebook# type: ignore
notebook = salome_notebook.NoteBook()

###
### GEOM component
###

import GEOM# type: ignore
from salome.geom import geomBuilder# type: ignore
import math
import SALOMEDS# type: ignore


geompy = geomBuilder.New()

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

upper_points = []
lower_points = []

with open(coordinates_file_name, 'r') as coordinates_file:
    airfoil_name = coordinates_file.readline()
    coordinates_number = coordinates_file.readline()
    upper_surface_coord_number = int(coordinates_number.split()[0])
    lower_surface_coord_number = int(coordinates_number.split()[1])
    print('Airfoil: ', airfoil_name)
    print('Number of upper points: ',  upper_surface_coord_number)
    print('Number of lower points: ', lower_surface_coord_number)
    vertex_id = 0
    lower_surface = False
    empty_line = coordinates_file.readline()
    line = coordinates_file.readline()
    while line:
        if '0.' in line:
            x = float(line.split()[0])
            y = float(line.split()[1])
            vertex = geompy.MakeVertex(x, y, 0)

            vertex_id += 1
            if lower_surface:
                lower_points.append(vertex)
            else:
                upper_points.append(vertex)
        if not line.strip():
            lower_surface = True
        line = coordinates_file.readline()

if len(upper_points) != upper_surface_coord_number or len(lower_points) != lower_surface_coord_number:
    print('Warning: number of points different than prescribed in input file:')
    print('Number of upper points (input): ',  upper_surface_coord_number)
    print('Number of upper points: ',  len(upper_points))
    print('Number of lower points (input): ', lower_surface_coord_number)
    print('Number of lower points: ', len(lower_points))

Upper_Airfoil = geompy.MakeInterpol(upper_points, False, False)
Lower_Airfoil = geompy.MakeInterpol(lower_points, False, False)

Upper_Airfoil_Divided = geompy.DivideEdge(Upper_Airfoil, -1, 0.5, 1)
Lower_Airfoil_Divided = geompy.DivideEdge(Lower_Airfoil, -1, 0.5, 1)

[Curve_UpperSurface_LE,Curve_UpperSurface_TE0] = geompy.ExtractShapes(Upper_Airfoil_Divided, geompy.ShapeType["EDGE"], True)
[Curve_LowerSurface_LE,Curve_LowerSurface_TE0] = geompy.ExtractShapes(Lower_Airfoil_Divided, geompy.ShapeType["EDGE"], True)

Curve_UpperSurface_TE = geompy.ChangeOrientationShellCopy(Curve_UpperSurface_TE0)
Curve_LowerSurface_TE = geompy.ChangeOrientationShellCopy(Curve_LowerSurface_TE0)

#Create face
Face_Airfoil_Before = geompy.MakeFaceWires([Curve_UpperSurface_LE, Curve_UpperSurface_TE, Curve_LowerSurface_TE, Curve_LowerSurface_LE], 1)
Face_Airfoil = geompy.MakeTranslation(Face_Airfoil_Before, -0.5, 0, 0)

#Rotate around center to AOA
geompy.Rotate(Face_Airfoil, OZ, -AOA*math.pi/180.0)

#Create domain
Face_Domain = geompy.MakeFaceHW(Domain_Length, Domain_Width, 1)

#Cut the airfoil from the domain
Cut_Domain = geompy.MakeCutList(Face_Domain, [Face_Airfoil], True)

#Explode edges
[Edge_1,Edge_2,Edge_LowerSurface_LE,Edge_UpperSurface_LE,Edge_LowerSurface_TE,Edge_UpperSurface_TE,Edge_7,Edge_8] = geompy.ExtractShapes(Cut_Domain, geompy.ShapeType["EDGE"], True)

[Auto_group_for_Sub_mesh_1,Auto_group_for_Sub_mesh_1_1] = geompy.ExtractShapes(Edge_LowerSurface_TE, geompy.ShapeType["VERTEX"], True)

#LowerSurface
Auto_group_for_Sub_mesh_1 = geompy.CreateGroup(Cut_Domain, geompy.ShapeType["EDGE"])
geompy.UnionList(Auto_group_for_Sub_mesh_1, [Edge_LowerSurface_LE, Edge_LowerSurface_TE])

#LowerSurface2
Auto_group_for_Sub_mesh_1_1 = geompy.CreateGroup(Cut_Domain, geompy.ShapeType["EDGE"])
geompy.UnionList(Auto_group_for_Sub_mesh_1_1, [Edge_LowerSurface_LE, Edge_LowerSurface_TE])

#UpperSurface
Auto_group_for_Sub_mesh_2 = geompy.CreateGroup(Cut_Domain, geompy.ShapeType["EDGE"])
geompy.UnionList(Auto_group_for_Sub_mesh_2, [Edge_UpperSurface_LE, Edge_UpperSurface_TE])

#Body
Body_Sub_mesh = geompy.CreateGroup(Cut_Domain, geompy.ShapeType["EDGE"])
geompy.UnionList(Body_Sub_mesh, [Edge_LowerSurface_LE, Edge_LowerSurface_TE, Edge_UpperSurface_LE, Edge_UpperSurface_TE])

#FarField
Auto_group_for_Sub_mesh_1_2 = geompy.CreateGroup(Cut_Domain, geompy.ShapeType["EDGE"])
geompy.UnionList(Auto_group_for_Sub_mesh_1_2, [Edge_1, Edge_2, Edge_7, Edge_8])

geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Upper_Airfoil, 'Upper_Airfoil' )
geompy.addToStudy( Lower_Airfoil, 'Lower_Airfoil' )
geompy.addToStudy( Upper_Airfoil_Divided, 'Upper_Airfoil_Divided' )
geompy.addToStudy( Lower_Airfoil_Divided, 'Lower_Airfoil_Divided' )
geompy.addToStudyInFather( Upper_Airfoil_Divided, Curve_UpperSurface_LE, 'Curve_UpperSurface_LE' )
geompy.addToStudyInFather( Upper_Airfoil_Divided, Curve_UpperSurface_TE, 'Curve_UpperSurface_TE' )
geompy.addToStudyInFather( Lower_Airfoil_Divided, Curve_LowerSurface_LE, 'Curve_LowerSurface_LE' )
geompy.addToStudyInFather( Lower_Airfoil_Divided, Curve_LowerSurface_TE, 'Curve_LowerSurface_TE' )
geompy.addToStudy( Face_Airfoil, 'Face_Airfoil' )
geompy.addToStudy( Face_Domain, 'Face_Domain' )
geompy.addToStudy( Cut_Domain, 'Cut_Domain' )
geompy.addToStudyInFather( Cut_Domain, Edge_LowerSurface_LE, 'Edge_LowerSurface_LE' )
geompy.addToStudyInFather( Cut_Domain, Edge_UpperSurface_LE, 'Edge_UpperSurface_LE' )
geompy.addToStudyInFather( Cut_Domain, Edge_LowerSurface_TE, 'Edge_LowerSurface_TE' )
geompy.addToStudyInFather( Cut_Domain, Edge_UpperSurface_TE, 'Edge_UpperSurface_TE' )
geompy.addToStudyInFather( Cut_Domain, Edge_7, 'Edge_7' )
geompy.addToStudyInFather( Cut_Domain, Edge_8, 'Edge_8' )
geompy.addToStudyInFather( Cut_Domain, Auto_group_for_Sub_mesh_1, 'Auto_group_for_Sub-mesh_1' )
geompy.addToStudyInFather( Cut_Domain, Auto_group_for_Sub_mesh_1_1, 'Auto_group_for_Sub-mesh_1' )
geompy.addToStudyInFather( Cut_Domain, Auto_group_for_Sub_mesh_2, 'Auto_group_for_Sub-mesh_2' )
geompy.addToStudyInFather( Cut_Domain, Auto_group_for_Sub_mesh_1_2, 'Auto_group_for_Sub-mesh_1' )

###
### SMESH component
###

import  SMESH, SALOMEDS# type: ignore
from salome.smesh import smeshBuilder # type: ignore

smesh = smeshBuilder.New()

#Set NETGEN
NETGEN_1D_2D = smesh.CreateHypothesis('NETGEN_2D', 'NETGENEngine')
NETGEN_2D_Parameters_1 = smesh.CreateHypothesis('NETGEN_Parameters_2D', 'NETGENEngine')
NETGEN_2D_Parameters_1.SetMaxSize( FarField_MeshSize )
NETGEN_2D_Parameters_1.SetSecondOrder( 0 )
NETGEN_2D_Parameters_1.SetOptimize( 1 )
NETGEN_2D_Parameters_1.SetFineness( 5 )
NETGEN_2D_Parameters_1.SetGrowthRate( Growth_Rate )
NETGEN_2D_Parameters_1.SetNbSegPerEdge( 3 )
NETGEN_2D_Parameters_1.SetNbSegPerRadius( 5 )
NETGEN_2D_Parameters_1.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_1.SetFuseEdges( 1 )
NETGEN_2D_Parameters_1.SetQuadAllowed( 0 )
Fluid = smesh.Mesh(Cut_Domain)
status = Fluid.AddHypothesis(NETGEN_2D_Parameters_1)
status = Fluid.AddHypothesis(NETGEN_1D_2D)
NETGEN_2D_Parameters_1.SetMinSize( 1e-15 )

#Set submeshes
#Body
Regular_1D = Fluid.Segment(geom=Body_Sub_mesh)
#Set geometric mesh
Geometric_Progression_1 = Regular_1D.GeometricProgression(Airfoil_MeshSize,Ratio,[])

#Set farfield mesh
Regular_1D_2 = Fluid.Segment(geom=Auto_group_for_Sub_mesh_1_2)
Local_Length_1 = Regular_1D_2.LocalLength(FarField_MeshSize,None,1e-07)

#Mesh
isDone = Fluid.Compute()
Body = Regular_1D.GetSubMesh()
FarField = Regular_1D_2.GetSubMesh()

NumberOfNodes = Fluid.NbNodes()
print('Information about surface mesh:')
print('Number of nodes       :', NumberOfNodes)

## Set names of Mesh objects
smesh.SetName(NETGEN_1D_2D, 'NETGEN 1D-2D')
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Geometric_Progression_1, 'Geometric Progression_1')
smesh.SetName(NETGEN_2D_Parameters_1, 'NETGEN 2D Parameters_1')
smesh.SetName(FarField, 'PotentialWallCondition2D_Far_field_Auto1')
smesh.SetName(Local_Length_1, 'Local Length_1')
smesh.SetName(Body, 'Body2D_Body')
smesh.SetName(Fluid.GetMesh(), 'Fluid')

# Saving file to open from salome's gui
file_name = salome_output_path + "/generate_cosine.hdf"
salome.myStudy.SaveAs(file_name, salome.myStudy, 0)


sys.path.append("../../../KratosSalomePlugin") # adding root folder of plugin to path
import create_kratos_input_tui# type: ignore

mesh_description_domain = { "elements"   : {"Triangle" : {"Element2D3N" : 0} } }
mesh_description_wall   = { "conditions" : {"Edge"     : {"WallCondition2D2N" : 0} } }

Fluid.ExportMED(f"{salome_output_path}/model_mesh_0.med")

meshes = [
        create_kratos_input_tui.SalomeMesh(Fluid, mesh_description_domain, "Parts_Parts_Auto1"),
        create_kratos_input_tui.SalomeMesh(FarField, mesh_description_wall, "PotentialWallCondition2D_Far_field_Auto1"),
        create_kratos_input_tui.SalomeMesh(Body, mesh_description_wall, "Body2D_Body"),
        ]

create_kratos_input_tui.CreateMdpaFile(meshes, f"{salome_output_path}/model_mesh_0")

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(True)