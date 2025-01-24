#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.10.0 with dump python functionality
###

import os
import sys
import salome # type: ignore

salome.salome_init()
study = salome.myStudy
import salome_notebook # type: ignore
notebook = salome_notebook.NoteBook()

WingAngle        = 0.0
########################
# DomainHeight     = 6.0 #X
# DomainWidth      = 3.0 #Y
# DomainLength     = 6.0 #Z
DomainHeight     = 25.0 #X
DomainWidth      = 25.0 #Y
DomainLength     = 25.0 #Z
########################
MeshGrowthRate   = 0.2
########################
MeshMaxSize      = 2.0
MeshMinSize      = 0.007
########################
MeshWingSize     = 0.009
MeshTEWingSize   = 0.007
MeshLEWingSize   = 0.007
MeshTipWingSize  = 0.007
MeshRootFoilSize = 0.009
########################
# Wake_angle         = 3.06
# Wake_length        = DomainWidth
# Wake_max_mesh_size = 1.0
# Wake_min_mesh_size = 0.01
# TE_Wing_mesh_size  = 0.01

# Dir                     = os.getcwd()
onera_geometry_igs_path = "onera_m6_geometry/Onera_M6_geometry.igs"
MeshOutPutName          = "Fluid.med"
# Wake_output_name        = "Wake.stl"

# if not os.path.exists(f'{Dir}/SalomeFiles'):
#   os.mkdir(f'{Dir}/SalomeFiles')


###
### GEOM component
###
from salome.geom import geomBuilder # type: ignore
import math


geompy = geomBuilder.New()

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
onera_geometry = geompy.ImportIGES(onera_geometry_igs_path)
geompy.Rotate(onera_geometry, OY, WingAngle*math.pi/180.0)
Face_1 = geompy.MakeFaceHW(DomainLength, DomainHeight, 3)
Extrusion_1 = geompy.MakePrismDXDYDZ(Face_1, 0, DomainWidth, 0)

Domain = geompy.MakeCutList(Extrusion_1, [onera_geometry], True)

Inlet = geompy.CreateGroup(Domain, geompy.ShapeType["FACE"])
geompy.UnionIDs(Inlet, [13])

Outlet = geompy.CreateGroup(Domain, geompy.ShapeType["FACE"])
geompy.UnionIDs(Outlet, [20])

Wing = geompy.CreateGroup(Domain, geompy.ShapeType["FACE"])
geompy.UnionIDs(Wing, [44, 75, 66, 57])

Walls = geompy.CreateGroup(Domain, geompy.ShapeType["FACE"])
geompy.UnionIDs(Walls, [3, 39, 42])

FarField = geompy.CreateGroup(Domain, geompy.ShapeType["FACE"])
geompy.UnionIDs(FarField, [42, 13, 3, 39, 20])

TE_Wing = geompy.CreateGroup(Domain, geompy.ShapeType["EDGE"])
geompy.UnionIDs(TE_Wing, [56])

LE_Wing = geompy.CreateGroup(Domain, geompy.ShapeType["EDGE"])
geompy.UnionIDs(LE_Wing, [46])

Tip_Wing = geompy.CreateGroup(Domain, geompy.ShapeType["EDGE"])
geompy.UnionIDs(Tip_Wing, [70, 74, 68, 72])

Root_Airfoil = geompy.CreateGroup(Domain, geompy.ShapeType["EDGE"])
geompy.UnionIDs(Root_Airfoil, [36, 38, 34, 31])

TE_Points = geompy.CreateGroup(Domain, geompy.ShapeType["VERTEX"])
geompy.UnionIDs(TE_Points, [33, 55])

[Inlet, Outlet, Wing, Walls, FarField, TE_Wing, LE_Wing, Tip_Wing, Root_Airfoil, TE_Points] = geompy.GetExistingSubObjects(Domain, False)

#### WAKE
# wake_base = geompy.MakePrismDXDYDZ(TE_Wing, Wake_length, 0, 0)
# aux_wake_base = geompy.MakePrismDXDYDZ(TE_Wing, -0.1, 0, 0)
# geompy.Rotate(wake_base, TE_Wing, Wake_angle*math.pi/180.0)
# geompy.Rotate(aux_wake_base, TE_Wing, -WingAngle*math.pi/180.0)
# wake_aux = geompy.MakeFuseList([wake_base, aux_wake_base], True, True)
# Left_Boundary = geompy.CreateGroup(wake_aux, geompy.ShapeType["EDGE"])
# geompy.UnionIDs(Left_Boundary, [7, 16])
# Right_Boundary = geompy.CreateGroup(wake_aux, geompy.ShapeType["EDGE"])
# geompy.UnionIDs(Right_Boundary, [4])
# Extrusion_Left = geompy.MakePrismDXDYDZ(Left_Boundary, 0, -0.1, 0)
# Extrusion_Right = geompy.MakePrismDXDYDZ(Right_Boundary, 0, 0.1, 0)
# wake = geompy.MakeFuseList([wake_aux, Extrusion_Left, Extrusion_Right], True, True)
##########################################################################

geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( onera_geometry, 'onera_geometry' )
########## WAKE
# geompy.addToStudy( wake_base, 'wake_base' )
# geompy.addToStudy( aux_wake_base, 'aux_wake_base' )
# geompy.addToStudy( wake_aux, 'wake_aux' )
# geompy.addToStudyInFather( wake_aux, Left_Boundary, 'Left_Boundary' )
# geompy.addToStudyInFather( wake_aux, Right_Boundary, 'Right_Boundary' )
# geompy.addToStudy( Extrusion_Left, 'Extrusion_Left' )
# geompy.addToStudy( Extrusion_Right, 'Extrusion_Right' )
# geompy.addToStudy( wake, 'wake' )
#########################################################################
geompy.addToStudy( Face_1, 'Face_1' )
geompy.addToStudy( Extrusion_1, 'Extrusion_1' )
geompy.addToStudy( Domain, 'Domain' )
geompy.addToStudyInFather( Domain, Inlet, 'Inlet' )
geompy.addToStudyInFather( Domain, Outlet, 'Outlet' )
geompy.addToStudyInFather( Domain, Wing, 'Wing' )
geompy.addToStudyInFather( Domain, Walls, 'Walls' )
geompy.addToStudyInFather( Domain, FarField, 'FarField' )
geompy.addToStudyInFather( Domain, TE_Wing, 'TE_Wing' )
geompy.addToStudyInFather( Domain, LE_Wing, 'LE_Wing' )
geompy.addToStudyInFather( Domain, Tip_Wing, 'Tip_Wing' )
geompy.addToStudyInFather( Domain, Root_Airfoil, 'Root_Airfoil' )
geompy.addToStudyInFather( Domain, TE_Points, 'TE_Points' )


###
### SMESH component
###

import  SMESH # type: ignore
from salome.smesh import smeshBuilder # type: ignore

sys.path.append("../../../../../KratosSalomePlugin") # adding root folder of plugin to path
import create_kratos_input_tui # type: ignore

mesh_description_domain  = { "elements"   : {"Tetra"    : {"Element3D4N"       : 0}}}
mesh_description_surface = { "conditions" : {"Triangle" : {"WallCondition3D3N" : 0}}}
mesh_description_wall    = { "conditions" : {"Edge"     : {"WallCondition2D2N" : 0}}}

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Fluid = smesh.Mesh(Domain,'Fluid')
NETGEN_1D_2D_3D = Fluid.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
NETGEN_3D_Parameters = NETGEN_1D_2D_3D.Parameters()
NETGEN_3D_Parameters.SetMaxSize( MeshMaxSize )
NETGEN_3D_Parameters.SetMinSize( MeshMinSize )
NETGEN_3D_Parameters.SetSecondOrder( 0 )
NETGEN_3D_Parameters.SetOptimize( 3 )
NETGEN_3D_Parameters.SetFineness( 5 )
NETGEN_3D_Parameters.SetGrowthRate( MeshGrowthRate )
NETGEN_3D_Parameters.SetNbSegPerEdge( 15 )
NETGEN_3D_Parameters.SetNbSegPerRadius( 2 )
NETGEN_3D_Parameters.SetChordalError( -1 )
NETGEN_3D_Parameters.SetChordalErrorEnabled( 1 )
NETGEN_3D_Parameters.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters.SetFuseEdges( 1 )
NETGEN_3D_Parameters.SetQuadAllowed( 0 )
NETGEN_3D_Parameters.SetLocalSizeOnShape(Wing, MeshWingSize)
NETGEN_3D_Parameters.SetLocalSizeOnShape(TE_Wing, MeshTEWingSize)
NETGEN_3D_Parameters.SetLocalSizeOnShape(LE_Wing, MeshLEWingSize)
NETGEN_3D_Parameters.SetLocalSizeOnShape(Tip_Wing, MeshTipWingSize)
NETGEN_3D_Parameters.SetLocalSizeOnShape(Root_Airfoil, MeshRootFoilSize)
Inlet_1 = Fluid.GroupOnGeom(Inlet,'Inlet',SMESH.FACE)
Outlet_1 = Fluid.GroupOnGeom(Outlet,'Outlet',SMESH.FACE)
Wing_1 = Fluid.GroupOnGeom(Wing,'Wing',SMESH.FACE)
Walls_1 = Fluid.GroupOnGeom(Walls,'Walls',SMESH.FACE)
FarField_1 = Fluid.GroupOnGeom(FarField,'FarField',SMESH.FACE)
TE_Wing_1 = Fluid.GroupOnGeom(TE_Wing,'TE_Wing',SMESH.EDGE)
LE_Wing_1 = Fluid.GroupOnGeom(LE_Wing,'LE_Wing',SMESH.EDGE)
Tip_Wing_1 = Fluid.GroupOnGeom(Tip_Wing,'Tip_Wing',SMESH.EDGE)
Root_Airfoil_1 = Fluid.GroupOnGeom(Root_Airfoil,'Root_Airfoil',SMESH.EDGE)
TE_Points_1 = Fluid.GroupOnGeom(TE_Points,'TE_Points',SMESH.NODE)
isDone = Fluid.Compute()

print(" Volume mesh READY")

[ Inlet_1, Outlet_1, Wing_1, Walls_1, FarField_1, TE_Wing_1, LE_Wing_1, Tip_Wing_1, Root_Airfoil_1, TE_Points_1] = Fluid.GetGroups()
Domain_1 = Fluid.GroupOnGeom(Domain,'Domain',SMESH.VOLUME)
Domain_2 = Fluid.GroupOnGeom(Domain,'Domain',SMESH.NODE)
smesh.SetName(Fluid, 'Fluid')

try:
  Fluid.ExportMED( MeshOutPutName, 0, 41, 1, Fluid, 1, [], '',-1, 1 )
  pass
except:
  print('ExportPartToMED() failed. Invalid file name?')


## Set names of Mesh objects
smesh.SetName(NETGEN_1D_2D_3D.GetAlgorithm(), 'NETGEN 1D-2D-3D')
smesh.SetName(NETGEN_3D_Parameters, 'NETGEN 3D Parameters')
smesh.SetName(Inlet_1, 'Inlet')
smesh.SetName(Outlet_1, 'Outlet')
smesh.SetName(Wing_1, 'Wing')
smesh.SetName(Walls_1, 'Walls')
smesh.SetName(FarField_1, 'FarField')
smesh.SetName(Fluid.GetMesh(), 'Fluid')
smesh.SetName(Domain_1, 'Domain')
smesh.SetName(TE_Wing_1, 'TE_Wing')
smesh.SetName(LE_Wing_1, 'LE_Wing')
smesh.SetName(Tip_Wing_1, 'Tip_Wing')
smesh.SetName(Root_Airfoil_1, 'Root_Airfoil')
smesh.SetName(TE_Points_1, 'TE_Points')
smesh.SetName(Domain_2, 'Domain')


# ########## WAKE
# Wake = smesh.Mesh(wake,'Wake')
# NETGEN_1D_2D = Wake.Triangle(algo=smeshBuilder.NETGEN_1D2D)
# NETGEN_2D_Parameters = NETGEN_1D_2D.Parameters()
# NETGEN_2D_Parameters.SetMaxSize( Wake_max_mesh_size )
# NETGEN_2D_Parameters.SetMinSize( Wake_min_mesh_size )
# NETGEN_2D_Parameters.SetSecondOrder( 0 )
# NETGEN_2D_Parameters.SetOptimize( 1 )
# NETGEN_2D_Parameters.SetFineness( 2 )
# NETGEN_2D_Parameters.SetChordalError( -1 )
# NETGEN_2D_Parameters.SetChordalErrorEnabled( 0 )
# NETGEN_2D_Parameters.SetUseSurfaceCurvature( 1 )
# NETGEN_2D_Parameters.SetFuseEdges( 1 )
# NETGEN_2D_Parameters.SetWorstElemMeasure( 0 )
# NETGEN_2D_Parameters.SetQuadAllowed( 0 )
# NETGEN_2D_Parameters.SetLocalSizeOnShape(TE_Wing, TE_Wing_mesh_size)
# NETGEN_2D_Parameters.SetUseDelauney( 128 )
# NETGEN_2D_Parameters.SetCheckChartBoundary( 3 )
# isDone = Wake.Compute()

# print(" Wake mesh READY")

# try:
#   Wake.ExportSTL( Wake_output_name, 1 )
#   pass
# except:
#   print('ExportSTL() failed. Invalid file name?')
# ## Set names of Mesh objects
# smesh.SetName(NETGEN_1D_2D.GetAlgorithm(), 'NETGEN 1D-2D')
# smesh.SetName(NETGEN_2D_Parameters, 'NETGEN 2D Parameters')
# smesh.SetName(Wake.GetMesh(), 'Wake')

#Save salome files
salome.myStudy.SaveAs("salome_model.hdf", study, False)


meshes = [
          create_kratos_input_tui.SalomeMesh(Fluid, mesh_description_domain, "Fluid"),
          create_kratos_input_tui.SalomeMesh(FarField_1, mesh_description_surface, "FarField"),
          create_kratos_input_tui.SalomeMesh(TE_Wing_1, mesh_description_surface, "TE_Wing"),
          create_kratos_input_tui.SalomeMesh(Wing_1, mesh_description_surface, "Wing"),
          ]

create_kratos_input_tui.CreateMdpaFile(meshes, f"Fluid")

##########################################################################

print(' Information about volume mesh:')
print(' Number of nodes       :', Fluid.NbNodes())
print(' Number of elements    :', Fluid.NbTetras())
# print(' Information about wake mesh:')
# print(' Number of nodes       :', Wake.NbNodes())
# print(' Number of elements    :', Wake.NbTriangles())

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
