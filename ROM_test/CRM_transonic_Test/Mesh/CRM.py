#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.11.0 with dump python functionality
###

import os
import sys
import salome# type: ignore

salome.salome_init()
study = salome.myStudy
import salome_notebook# type: ignore
notebook = salome_notebook.NoteBook()

PlaneAngle = 0.0
########################
DomainHeight = 250.0 #X
DomainWidth  = 125.0 #Y
DomainLength = 250.0 #Z
########################
MeshGrowthRate = 0.2
########################
MeshMaxSize = 50.0
MeshMinSize = 0.05
########################
MeshPlaneMaxSize = 0.05
MeshPlaneMinSize = 0.05
########################
Wake_angle         = 2.31
Wake_length        = DomainWidth
Wake_max_mesh_size = 50.0
Wake_min_mesh_size = 0.1
TE_Plane_mesh_size = 0.1

Dir                   = os.getcwd()
crm_geometry_igs_path = Dir + "/crm_geometry/nasa_crm_original_TE_00.igs"
MeshOutPutName        = Dir + "/Fluid.med"
Wake_output_name      = Dir + "/wake_sup.stl"


###
### GEOM component
###

import GEOM # type: ignore
from salome.geom import geomBuilder# type: ignore
import math
import SALOMEDS # type: ignore

print('Init geometry load ...')
geompy = geomBuilder.New()

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
SOLID_1_ = geompy.ImportIGES(crm_geometry_igs_path, True)
print('finish geometry load')
print('Init scale ...')
CRM = geompy.MakeScaleTransform(SOLID_1_, None, 0.025)
print('Finished scale')
geompy.Rotate(CRM, OY, PlaneAngle*math.pi/180.0)
Face_1 = geompy.MakeFaceHW(DomainLength, DomainHeight, 3)
Extrusion_1 = geompy.MakePrismDXDYDZ(Face_1, 0, DomainWidth, 0)
Domain = geompy.MakeCutList(Extrusion_1, [CRM])
Inlet = geompy.CreateGroup(Domain, geompy.ShapeType["FACE"])
geompy.UnionIDs(Inlet, [13])
Walls = geompy.CreateGroup(Domain, geompy.ShapeType["FACE"])
geompy.UnionIDs(Walls, [74, 3, 71])
FarField = geompy.CreateGroup(Domain, geompy.ShapeType["FACE"])
geompy.UnionIDs(FarField, [71, 74, 3, 20, 13])
Plane = geompy.CreateGroup(Domain, geompy.ShapeType["FACE"])
geompy.UnionIDs(Plane, [228, 128, 167, 268, 175, 223, 185, 170, 123, 81, 278, 164, 249, 284, 90, 187, 263, 233, 172, 281, 215, 258, 108, 154, 220, 95, 133, 273, 159, 203, 76, 236, 208, 111, 238, 196, 118, 145, 102])
Plane_TE_Sup = geompy.CreateGroup(Domain, geompy.ShapeType["EDGE"])
geompy.UnionIDs(Plane_TE_Sup, [219, 198, 257])

Wing_Root_Point_Sup = geompy.CreateGroup(Domain, geompy.ShapeType["VERTEX"])
geompy.UnionIDs(Wing_Root_Point_Sup, [183])
Tail_Root_Point_Sup = geompy.CreateGroup(Domain, geompy.ShapeType["VERTEX"])
geompy.UnionIDs(Tail_Root_Point_Sup, [143])

[Inlet, Walls, FarField, Plane, Plane_TE_Sup, Wing_Root_Point_Sup, Tail_Root_Point_Sup] = geompy.GetExistingSubObjects(Domain, False)

Wing_Root_Line_Sup = geompy.MakePrismDXDYDZ(Wing_Root_Point_Sup, 0, -4, 0)
Wing_Root_Tail_Sup = geompy.MakePrismDXDYDZ(Tail_Root_Point_Sup, 0, -2, 0)
Plane_TE_Lines = geompy.MakeFuseList([Wing_Root_Line_Sup, Wing_Root_Tail_Sup, Plane_TE_Sup], True, True)
Wake_sup = geompy.MakePrismDXDYDZ(Plane_TE_Lines, DomainWidth, 0, DomainWidth*math.tan(Wake_angle*math.pi/180))

geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( SOLID_1_, 'SOLID(1)' )
geompy.addToStudy( CRM, 'CRM' )
geompy.addToStudy( Face_1, 'Face_1' )
geompy.addToStudy( Extrusion_1, 'Extrusion_1' )
geompy.addToStudy( Domain, 'Domain' )
geompy.addToStudyInFather( Domain, Inlet, 'Inlet' )
geompy.addToStudyInFather( Domain, Walls, 'Walls' )
geompy.addToStudyInFather( Domain, FarField, 'FarField' )
geompy.addToStudyInFather( Domain, Plane, 'Plane' )
geompy.addToStudyInFather( Domain, Plane_TE_Sup, 'Plane_TE_Sup' )
geompy.addToStudyInFather( Domain, Wing_Root_Point_Sup, 'Wing_Root_Point_Sup' )
geompy.addToStudyInFather( Domain, Tail_Root_Point_Sup, 'Tail_Root_Point_Sup' )
geompy.addToStudy( Wing_Root_Line_Sup, 'Wing_Root_Line_Sup' )
geompy.addToStudy( Wing_Root_Tail_Sup, 'Wing_Root_Tail_Sup' )
geompy.addToStudy( Plane_TE_Lines, 'Plane_TE_Lines' )
geompy.addToStudy( Wake_sup, 'Wake_sup' )

print('Finished geometry creation')
###
### SMESH component
###

import  SMESH, SALOMEDS # type: ignore
from salome.smesh import smeshBuilder # type: ignore

sys.path.append("../../../../KratosSalomePlugin") # adding root folder of plugin to path
import create_kratos_input_tui # type: ignore

mesh_description_domain  = { "elements"   : {"Tetra"    : {"Element3D4N"       : 0}}}
mesh_description_surface = { "conditions" : {"Triangle" : {"WallCondition3D3N" : 0}}}
mesh_description_wall    = { "conditions" : {"Edge"     : {"WallCondition2D2N" : 0}}}

print('Init meshing')

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Fluid = smesh.Mesh(Domain,'Fluid')
NETGEN_1D_2D_3D = Fluid.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
NETGEN_3D_Parameters = NETGEN_1D_2D_3D.Parameters()
NETGEN_3D_Parameters.SetMaxSize( MeshMaxSize )
NETGEN_3D_Parameters.SetMinSize( MeshMinSize )
NETGEN_3D_Parameters.SetSecondOrder( 0 )
NETGEN_3D_Parameters.SetOptimize( 1 )
NETGEN_3D_Parameters.SetFineness( 0 )
NETGEN_3D_Parameters.SetGrowthRate( MeshGrowthRate )
NETGEN_3D_Parameters.SetChordalError( -1 )
NETGEN_3D_Parameters.SetChordalErrorEnabled( 0 )
NETGEN_3D_Parameters.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters.SetFuseEdges( 1 )
NETGEN_3D_Parameters.SetQuadAllowed( 0 )
NETGEN_3D_Parameters.SetCheckChartBoundary( 3 )

Inlet_1 = Fluid.GroupOnGeom(Inlet,'Inlet',SMESH.FACE)
Walls_1 = Fluid.GroupOnGeom(Walls,'Walls',SMESH.FACE)
FarField_1 = Fluid.GroupOnGeom(FarField,'FarField',SMESH.FACE)
Plane_1 = Fluid.GroupOnGeom(Plane,'Plane',SMESH.FACE)
Plane_TE_Sup_1 = Fluid.GroupOnGeom(Plane_TE_Sup,'Plane_TE_Sup',SMESH.EDGE)
Wing_Root_Point_Sup_1 = Fluid.GroupOnGeom(Wing_Root_Point_Sup,'Wing_Root_Point_Sup',SMESH.NODE)
Tail_Root_Point_Sup_1 = Fluid.GroupOnGeom(Tail_Root_Point_Sup,'Tail_Root_Point_Sup',SMESH.NODE)

NETGEN_1D_2D = Fluid.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=Plane)
NETGEN_2D_Parameters_Plane = NETGEN_1D_2D.Parameters()
NETGEN_2D_Parameters_Plane.SetMaxSize( MeshPlaneMaxSize )
NETGEN_2D_Parameters_Plane.SetMinSize( MeshPlaneMinSize )
NETGEN_2D_Parameters_Plane.SetSecondOrder( 0 )
NETGEN_2D_Parameters_Plane.SetOptimize( 1 )
NETGEN_2D_Parameters_Plane.SetFineness( 1 )
NETGEN_2D_Parameters_Plane.SetChordalError( -1 )
NETGEN_2D_Parameters_Plane.SetChordalErrorEnabled( 0 )
NETGEN_2D_Parameters_Plane.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_Plane.SetFuseEdges( 1 )
NETGEN_2D_Parameters_Plane.SetQuadAllowed( 0 )
NETGEN_2D_Parameters_Plane.SetWorstElemMeasure( 28911 )
NETGEN_2D_Parameters_Plane.SetUseDelauney( 96 )
NETGEN_2D_Parameters_Plane.SetCheckChartBoundary( 3 )

Wake = smesh.Mesh(Wake_sup,'Wake')
NETGEN_1D_2D_1 = Wake.Triangle(algo=smeshBuilder.NETGEN_1D2D)
NETGEN_2D_Parameters_Wake = NETGEN_1D_2D_1.Parameters()
NETGEN_2D_Parameters_Wake.SetMaxSize( MeshMaxSize )
NETGEN_2D_Parameters_Wake.SetMinSize( MeshMinSize )
NETGEN_2D_Parameters_Wake.SetSecondOrder( 0 )
NETGEN_2D_Parameters_Wake.SetOptimize( 1 )
NETGEN_2D_Parameters_Wake.SetFineness( 2 )
NETGEN_2D_Parameters_Wake.SetChordalError( -1 )
NETGEN_2D_Parameters_Wake.SetChordalErrorEnabled( 0 )
NETGEN_2D_Parameters_Wake.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_Wake.SetFuseEdges( 1 )
NETGEN_2D_Parameters_Wake.SetQuadAllowed( 0 )
NETGEN_2D_Parameters_Wake.SetLocalSizeOnShape(Plane_TE_Lines, MeshMinSize)
NETGEN_2D_Parameters_Wake.SetWorstElemMeasure( 24392 )
NETGEN_2D_Parameters_Wake.SetUseDelauney( 96 )
NETGEN_2D_Parameters_Wake.SetCheckChartBoundary( 3 )

isDone = Fluid.Compute()

[ Inlet_1, Walls_1, FarField_1, Plane_1, Plane_TE_Sup_1, Wing_Root_Point_Sup_1, Tail_Root_Point_Sup_1 ] = Fluid.GetGroups()

isDone = Wake.Compute()

Domain_1 = Fluid.GroupOnGeom(Domain,'Domain',SMESH.VOLUME)
smesh.SetName(Fluid, 'Fluid')

try:
  Fluid.ExportMED( MeshOutPutName, 0, 41, 1, Fluid, 1, [], '',-1, 1 )
  pass
except:
  print('ExportPartToMED() failed. Invalid file name?')

try:
  Wake.ExportSTL( Wake_output_name, 1 )
  pass
except:
  print('ExportSTL() failed. Invalid file name?')
Sub_mesh_Plane = NETGEN_1D_2D.GetSubMesh()


## Set names of Mesh objects
smesh.SetName(Tail_Root_Point_Sup_1, 'Tail_Root_Point_Sup')
smesh.SetName(Wing_Root_Point_Sup_1, 'Wing_Root_Point_Sup')
smesh.SetName(Sub_mesh_Plane, 'Sub-mesh_Plane')
smesh.SetName(NETGEN_1D_2D_3D.GetAlgorithm(), 'NETGEN 1D-2D-3D')
smesh.SetName(NETGEN_1D_2D.GetAlgorithm(), 'NETGEN 1D-2D')
smesh.SetName(Inlet_1, 'Inlet')
smesh.SetName(Domain_1, 'Domain')
smesh.SetName(FarField_1, 'FarField')
smesh.SetName(Walls_1, 'Walls')
smesh.SetName(Plane_1, 'Plane')
smesh.SetName(Fluid.GetMesh(), 'Fluid')
smesh.SetName(Wake.GetMesh(), 'Wake')
smesh.SetName(Plane_TE_Sup_1, 'Plane_TE_Sup')
smesh.SetName(NETGEN_2D_Parameters_Wake, 'NETGEN 2D Parameters_Wake')
smesh.SetName(NETGEN_2D_Parameters_Plane, 'NETGEN 2D Parameters_Plane')
smesh.SetName(NETGEN_3D_Parameters, 'NETGEN 3D Parameters')

#Save salome files
salome.myStudy.SaveAs(Dir + "/salome_model.hdf", study, False)

meshes = [
          create_kratos_input_tui.SalomeMesh(Fluid, mesh_description_domain, "Domain"),
          create_kratos_input_tui.SalomeMesh(FarField_1, mesh_description_surface, "FarField"),
          create_kratos_input_tui.SalomeMesh(Plane_TE_Sup_1, mesh_description_surface, "Plane_TE"),
          create_kratos_input_tui.SalomeMesh(Plane_1, mesh_description_surface, "Plane"),
          ]

create_kratos_input_tui.CreateMdpaFile(meshes, f"{Dir}/Fluid")

##########################################################################

print(' Information about volume mesh:')
print(' Number of nodes       :', Fluid.NbNodes())
print(' Number of elements    :', Fluid.NbTetras())
print(' Information about wake mesh:')
print(' Number of nodes       :', Wake.NbNodes())
print(' Number of elements    :', Wake.NbTriangles())

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
