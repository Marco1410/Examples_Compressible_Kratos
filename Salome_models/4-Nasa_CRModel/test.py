#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.11.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/media/marco-kratos-pc/datos/Examples_Compressible_Kratos/Salome_models/nasa_model_mod')

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS


geompy = geomBuilder.New()

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
nasa_crm_igs_1 = geompy.ImportIGES("/home/marco-kratos-pc/Descargas/nasa_crm.igs")
[SOLID_1_] = geompy.ExtractShapes(nasa_crm_igs_1, geompy.ShapeType["SOLID"], True)
Face_1 = geompy.MakeFaceHW(40, 40, 3)
Extrusion_1 = geompy.MakePrismDXDYDZ(Face_1, 0, 20, 0)
geompy.Rotate(nasa_crm_igs_1, OY, 2.31*math.pi/180.0)
Domain = geompy.MakeCutList(Extrusion_1, [SOLID_1_])
Inlet = geompy.CreateGroup(Domain, geompy.ShapeType["FACE"])
geompy.UnionIDs(Inlet, [13])
Outlet = geompy.CreateGroup(Domain, geompy.ShapeType["FACE"])
geompy.UnionIDs(Outlet, [20])
Walls = geompy.CreateGroup(Domain, geompy.ShapeType["FACE"])
geompy.UnionIDs(Walls, [60, 57, 3, 27])
FarField = geompy.CreateGroup(Domain, geompy.ShapeType["FACE"])
geompy.UnionIDs(FarField, [27, 60, 3, 57, 13, 20])
Plane = geompy.CreateGroup(Domain, geompy.ShapeType["FACE"])
geompy.UnionIDs(Plane, [127, 179, 124, 146, 204, 189, 221, 95, 166, 224, 108, 199, 159, 62, 174, 113, 202, 88, 211, 102, 79, 141, 149, 214, 169, 229, 132, 184, 194, 156, 74])
Wake_TE = geompy.CreateGroup(Domain, geompy.ShapeType["EDGE"])
geompy.UnionIDs(Wake_TE, [220, 136, 210])
[Inlet, Outlet, Walls, FarField, Plane, Wake_TE] = geompy.GetExistingSubObjects(Domain, False)
wake_lines = geompy.CreateGroup(SOLID_1_, geompy.ShapeType["EDGE"])
geompy.UnionIDs(wake_lines, [189, 179, 25])
Wake = geompy.MakePrismDXDYDZ(wake_lines, 20, 0, 0)
geompy.DifferenceIDs(FarField, [27, 60, 3, 57, 13, 20])
geompy.UnionIDs(FarField, [60, 3, 57, 13, 20])
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( nasa_crm_igs_1, 'nasa_crm.igs_1' )
geompy.addToStudy( Face_1, 'Face_1' )
geompy.addToStudy( Extrusion_1, 'Extrusion_1' )
geompy.addToStudyInFather( nasa_crm_igs_1, SOLID_1_, 'SOLID(1)' )
geompy.addToStudy( Domain, 'Domain' )
geompy.addToStudyInFather( Domain, Inlet, 'Inlet' )
geompy.addToStudyInFather( Domain, Outlet, 'Outlet' )
geompy.addToStudyInFather( Domain, Walls, 'Walls' )
geompy.addToStudyInFather( Domain, FarField, 'FarField' )
geompy.addToStudyInFather( Domain, Plane, 'Plane' )
geompy.addToStudyInFather( Domain, Wake_TE, 'Wake_TE' )
geompy.addToStudyInFather( SOLID_1_, wake_lines, 'wake_lines' )
geompy.addToStudy( Wake, 'Wake' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Fluid = smesh.Mesh(Domain,'Fluid')
NETGEN_1D_2D_3D = Fluid.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
NETGEN_3D_Parameters = NETGEN_1D_2D_3D.Parameters()
Inlet_1 = Fluid.GroupOnGeom(Inlet,'Inlet',SMESH.FACE)
Outlet_1 = Fluid.GroupOnGeom(Outlet,'Outlet',SMESH.FACE)
Walls_1 = Fluid.GroupOnGeom(Walls,'Walls',SMESH.FACE)
FarField_1 = Fluid.GroupOnGeom(FarField,'FarField',SMESH.FACE)
Plane_1 = Fluid.GroupOnGeom(Plane,'Plane',SMESH.FACE)
Wake_TE_1 = Fluid.GroupOnGeom(Wake_TE,'Wake_TE',SMESH.EDGE)
Regular_1D = Fluid.Segment(geom=Plane)
Local_Length_Plane = Regular_1D.LocalLength(0.06,None,1e-07)
NETGEN_2D = Fluid.Triangle(algo=smeshBuilder.NETGEN_2D,geom=Plane)
NETGEN_2D_Parameters_Plane = NETGEN_2D.Parameters()
NETGEN_2D_Parameters_Plane.SetWorstElemMeasure( 0 )
[ Inlet_1, Outlet_1, Walls_1, FarField_1, Plane_1, Wake_TE_1 ] = Fluid.GetGroups()
Wake_1 = smesh.Mesh(Wake,'Wake')
NETGEN_1D_2D = Wake_1.Triangle(algo=smeshBuilder.NETGEN_1D2D)
NETGEN_2D_Parameters_Wake = NETGEN_1D_2D.Parameters()
smesh.SetName(Fluid, 'Fluid')
try:
  Fluid.ExportMED( r'/home/marco-kratos-pc/Descargas/Mesh_Domain.med', 0, 41, 1, Fluid, 1, [], '',-1, 1 )
  pass
except:
  print('ExportPartToMED() failed. Invalid file name?')
try:
  Wake_1.ExportSTL( r'/home/marco-kratos-pc/Descargas/Wake.stl', 1 )
  pass
except:
  print('ExportSTL() failed. Invalid file name?')
Domain_1 = Fluid.GroupOnGeom(Domain,'Domain',SMESH.VOLUME)
Domain_2 = Fluid.GroupOnGeom(Domain,'Domain',SMESH.NODE)
smesh.SetName(Fluid, 'Fluid')
try:
  Fluid.ExportMED( r'/home/marco-kratos-pc/Descargas/Mesh_Domain.med', 0, 41, 1, Fluid, 1, [], '',-1, 1 )
  pass
except:
  print('ExportPartToMED() failed. Invalid file name?')
[ Inlet_1, Outlet_1, Walls_1, FarField_1, Plane_1, Wake_TE_1, Domain_1, Domain_2 ] = Fluid.GetGroups()
NETGEN_2D_Parameters_Wake.SetLocalSizeOnShape(wake_lines, 0.001)
smesh.SetName(Fluid, 'Fluid')
try:
  Fluid.ExportMED( r'/home/marco-kratos-pc/Descargas/Mesh_Domain.med', 0, 41, 1, Fluid, 1, [], '',-1, 1 )
  pass
except:
  print('ExportPartToMED() failed. Invalid file name?')
try:
  Wake_1.ExportSTL( r'/home/marco-kratos-pc/Descargas/Wake.stl', 1 )
  pass
except:
  print('ExportSTL() failed. Invalid file name?')
[ Inlet_1, Outlet_1, Walls_1, FarField_1, Plane_1, Wake_TE_1, Domain_1, Domain_2 ] = Fluid.GetGroups()
smesh.SetName(Fluid, 'Fluid')
try:
  Fluid.ExportMED( r'/home/marco-kratos-pc/Descargas/Mesh_Domain.med', 0, 41, 1, Fluid, 1, [], '',-1, 1 )
  pass
except:
  print('ExportPartToMED() failed. Invalid file name?')
[ Inlet_1, Outlet_1, Walls_1, FarField_1, Plane_1, Wake_TE_1, Domain_1, Domain_2 ] = Fluid.GetGroups()
smesh.SetName(Fluid, 'Fluid')
try:
  Fluid.ExportMED( r'/home/marco-kratos-pc/Descargas/Mesh_Domain.med', 0, 41, 1, Fluid, 1, [], '',-1, 1 )
  pass
except:
  print('ExportPartToMED() failed. Invalid file name?')
[ Inlet_1, Outlet_1, Walls_1, FarField_1, Plane_1, Wake_TE_1, Domain_1, Domain_2 ] = Fluid.GetGroups()
smesh.SetName(Fluid, 'Fluid')
try:
  Fluid.ExportMED( r'/home/marco-kratos-pc/Descargas/Mesh_Domain.med', 0, 41, 1, Fluid, 1, [], '',-1, 1 )
  pass
except:
  print('ExportPartToMED() failed. Invalid file name?')
[ Inlet_1, Outlet_1, Walls_1, FarField_1, Plane_1, Wake_TE_1, Domain_1, Domain_2 ] = Fluid.GetGroups()
NETGEN_3D_Parameters.SetMaxSize( 0.5 )
NETGEN_3D_Parameters.SetSecondOrder( 0 )
NETGEN_3D_Parameters.SetOptimize( 1 )
NETGEN_3D_Parameters.SetFineness( 4 )
NETGEN_3D_Parameters.SetChordalError( 0 )
NETGEN_3D_Parameters.SetChordalErrorEnabled( 0 )
NETGEN_3D_Parameters.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters.SetFuseEdges( 1 )
NETGEN_3D_Parameters.SetQuadAllowed( 0 )
NETGEN_2D_Parameters_Plane.SetOptimize( 1 )
NETGEN_2D_Parameters_Plane.SetFineness( 4 )
NETGEN_2D_Parameters_Plane.SetChordalError( 0 )
NETGEN_2D_Parameters_Plane.SetChordalErrorEnabled( 0 )
NETGEN_2D_Parameters_Plane.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_Plane.SetQuadAllowed( 0 )
NETGEN_2D_Parameters_Wake.SetMaxSize( 0.5 )
NETGEN_2D_Parameters_Wake.SetSecondOrder( 0 )
NETGEN_2D_Parameters_Wake.SetOptimize( 1 )
NETGEN_2D_Parameters_Wake.SetFineness( 4 )
NETGEN_2D_Parameters_Wake.SetChordalError( 0 )
NETGEN_2D_Parameters_Wake.SetChordalErrorEnabled( 0 )
NETGEN_2D_Parameters_Wake.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_Wake.SetFuseEdges( 1 )
NETGEN_2D_Parameters_Wake.SetWorstElemMeasure( 0 )
NETGEN_2D_Parameters_Wake.SetQuadAllowed( 0 )
NETGEN_2D_Parameters_Wake.SetLocalSizeOnShape(wake_lines, 0.001)
[ Inlet_1, Outlet_1, Walls_1, FarField_1, Plane_1, Wake_TE_1, Domain_1, Domain_2 ] = Fluid.GetGroups()
smesh.SetName(Fluid, 'Fluid')
try:
  Fluid.ExportMED( r'/home/marco-kratos-pc/Descargas/Mesh_Domain.med', 0, 41, 1, Fluid, 1, [], '',-1, 1 )
  pass
except:
  print('ExportPartToMED() failed. Invalid file name?')
try:
  Wake_1.ExportSTL( r'/home/marco-kratos-pc/Descargas/Wake.stl', 1 )
  pass
except:
  print('ExportSTL() failed. Invalid file name?')
[ Inlet_1, Outlet_1, Walls_1, FarField_1, Plane_1, Wake_TE_1, Domain_1, Domain_2 ] = Fluid.GetGroups()
smesh.SetName(Fluid, 'Fluid')
try:
  Fluid.ExportMED( r'/home/marco-kratos-pc/Descargas/Mesh_Domain.med', 0, 41, 1, Fluid, 1, [], '',-1, 1 )
  pass
except:
  print('ExportPartToMED() failed. Invalid file name?')
NETGEN_2D_Parameters_Wake.SetMinSize( 0.001 )
NETGEN_2D_Parameters_Wake.SetUseDelauney( 128 )
NETGEN_2D_Parameters_Wake.SetCheckChartBoundary( 3 )
NETGEN_2D_Parameters_Plane.SetMaxSize( 0.01 )
[ Inlet_1, Outlet_1, Walls_1, FarField_1, Plane_1, Wake_TE_1, Domain_1, Domain_2 ] = Fluid.GetGroups()
NETGEN_3D_Parameters.SetMinSize( 0.0025 )
NETGEN_3D_Parameters.SetCheckChartBoundary( 3 )
[ Inlet_1, Outlet_1, Walls_1, FarField_1, Plane_1, Wake_TE_1, Domain_1, Domain_2 ] = Fluid.GetGroups()
NETGEN_2D_Parameters_Plane.SetMinSize( 0.01 )
NETGEN_2D_Parameters_Plane.SetUseDelauney( 128 )
NETGEN_2D_Parameters_Plane.SetCheckChartBoundary( 3 )
Local_Length_Plane.SetLength( 0.01 )
Local_Length_Plane.SetPrecision( 1e-07 )
isDone = Fluid.Compute()
[ Inlet_1, Outlet_1, Walls_1, FarField_1, Plane_1, Wake_TE_1, Domain_1, Domain_2 ] = Fluid.GetGroups()
smesh.SetName(Fluid, 'Fluid')
try:
  Fluid.ExportMED( r'/home/marco-kratos-pc/Descargas/Mesh_Domain.med', 0, 41, 1, Fluid, 1, [], '',-1, 1 )
  pass
except:
  print('ExportPartToMED() failed. Invalid file name?')
[ Inlet_1, Outlet_1, Walls_1, FarField_1, Plane_1, Wake_TE_1, Domain_1, Domain_2 ] = Fluid.GetGroups()
smesh.SetName(Fluid, 'Fluid')
try:
  Fluid.ExportMED( r'/media/marco-kratos-pc/datos/Examples_Compressible_Kratos/Salome_models/nasa_model_mod/Mesh_Domain.med', 0, 41, 1, Fluid, 1, [], '',-1, 1 )
  pass
except:
  print('ExportPartToMED() failed. Invalid file name?')
isDone = Wake_1.Compute()
Sub_mesh_Plane = Regular_1D.GetSubMesh()


## Set names of Mesh objects
smesh.SetName(NETGEN_1D_2D_3D.GetAlgorithm(), 'NETGEN 1D-2D-3D')
smesh.SetName(NETGEN_2D.GetAlgorithm(), 'NETGEN 2D')
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(NETGEN_2D_Parameters_Plane, 'NETGEN 2D Parameters_Plane')
smesh.SetName(NETGEN_1D_2D.GetAlgorithm(), 'NETGEN 1D-2D')
smesh.SetName(Local_Length_Plane, 'Local Length_Plane')
smesh.SetName(NETGEN_3D_Parameters, 'NETGEN 3D Parameters')
smesh.SetName(NETGEN_2D_Parameters_Wake, 'NETGEN 2D Parameters_Wake')
smesh.SetName(Inlet_1, 'Inlet')
smesh.SetName(Outlet_1, 'Outlet')
smesh.SetName(Walls_1, 'Walls')
smesh.SetName(FarField_1, 'FarField')
smesh.SetName(Plane_1, 'Plane')
smesh.SetName(Fluid.GetMesh(), 'Fluid')
smesh.SetName(Wake_1.GetMesh(), 'Wake')
smesh.SetName(Sub_mesh_Plane, 'Sub-mesh_Plane')
smesh.SetName(Domain_1, 'Domain')
smesh.SetName(Wake_TE_1, 'Wake_TE')
smesh.SetName(Domain_2, 'Domain')

###
### PARAVIS component
###


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
