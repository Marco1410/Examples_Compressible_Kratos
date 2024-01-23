#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.11.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/media/marco-kratos-pc/datos/Examples_Compressible_Kratos/Salome_models/3-Onera/onera_m6_geometry')

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
SOLID_1_ = geompy.ImportIGES("/media/marco-kratos-pc/datos/Examples_Compressible_Kratos/Salome_models/3-Onera/onera_m6_geometry/Onera_M6_geometry.igs")
Wake_TE = geompy.CreateGroup(SOLID_1_, geompy.ShapeType["EDGE"])
geompy.UnionIDs(Wake_TE, [16])
Wake = geompy.MakePrismDXDYDZ(Wake_TE, 15, 0, 0)
Face_1 = geompy.MakeFaceHW(25, 25, 3)
Extrusion = geompy.MakePrismDXDYDZ(Face_1, 0, 12, 0)
Domain = geompy.MakeCutList(Extrusion, [SOLID_1_])
Inlet = geompy.CreateGroup(Domain, geompy.ShapeType["FACE"])
geompy.UnionIDs(Inlet, [13])
Outlet = geompy.CreateGroup(Domain, geompy.ShapeType["FACE"])
geompy.UnionIDs(Outlet, [20])
FarField = geompy.CreateGroup(Domain, geompy.ShapeType["FACE"])
geompy.UnionIDs(FarField, [13, 20, 42, 39, 3])
Walls = geompy.CreateGroup(Domain, geompy.ShapeType["FACE"])
geompy.UnionIDs(Walls, [3, 39, 42])
Wake_TE_1 = geompy.CreateGroup(Domain, geompy.ShapeType["EDGE"])
geompy.UnionIDs(Wake_TE_1, [56])
Wing = geompy.CreateGroup(Domain, geompy.ShapeType["FACE"])
geompy.UnionIDs(Wing, [75, 66, 44, 57])
[Inlet, Outlet, FarField, Walls, Wake_TE_1, Wing] = geompy.GetExistingSubObjects(Domain, False)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( SOLID_1_, 'SOLID(1)' )
geompy.addToStudyInFather( SOLID_1_, Wake_TE, 'Wake_TE' )
geompy.addToStudy( Wake, 'Wake' )
geompy.addToStudy( Face_1, 'Face_1' )
geompy.addToStudy( Extrusion, 'Extrusion' )
geompy.addToStudy( Domain, 'Domain' )
geompy.addToStudyInFather( Domain, Inlet, 'Inlet' )
geompy.addToStudyInFather( Domain, Outlet, 'Outlet' )
geompy.addToStudyInFather( Domain, FarField, 'FarField' )
geompy.addToStudyInFather( Domain, Walls, 'Walls' )
geompy.addToStudyInFather( Domain, Wake_TE_1, 'Wake_TE' )
geompy.addToStudyInFather( Domain, Wing, 'Wing' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()

Fluid = smesh.Mesh(Domain,'Fluid')
NETGEN_1D_2D_3D = Fluid.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
NETGEN_3D_Parameters_1 = NETGEN_1D_2D_3D.Parameters()
NETGEN_3D_Parameters_1.SetMaxSize( 0.5 )
NETGEN_3D_Parameters_1.SetMinSize( 0.04 )
NETGEN_3D_Parameters_1.SetSecondOrder( 0 )
NETGEN_3D_Parameters_1.SetOptimize( 1 )
NETGEN_3D_Parameters_1.SetFineness( 4 )
NETGEN_3D_Parameters_1.SetChordalError( -1 )
NETGEN_3D_Parameters_1.SetChordalErrorEnabled( 0 )
NETGEN_3D_Parameters_1.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_1.SetFuseEdges( 1 )
NETGEN_3D_Parameters_1.SetQuadAllowed( 0 )
NETGEN_3D_Parameters_1.SetCheckChartBoundary( 3 )
Inlet_1 = Fluid.GroupOnGeom(Inlet,'Inlet',SMESH.FACE)
Outlet_1 = Fluid.GroupOnGeom(Outlet,'Outlet',SMESH.FACE)
FarField_1 = Fluid.GroupOnGeom(FarField,'FarField',SMESH.FACE)
Walls_1 = Fluid.GroupOnGeom(Walls,'Walls',SMESH.FACE)
Wake_TE_2 = Fluid.GroupOnGeom(Wake_TE_1,'Wake_TE',SMESH.EDGE)
Wing_1 = Fluid.GroupOnGeom(Wing,'Wing',SMESH.FACE)
NETGEN_1D_2D = Fluid.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=Wing)
NETGEN_2D_Parameters_Wing = NETGEN_1D_2D.Parameters()
NETGEN_2D_Parameters_Wing.SetMaxSize( 0.04 )
NETGEN_2D_Parameters_Wing.SetMinSize( 0.04 )
NETGEN_2D_Parameters_Wing.SetSecondOrder( 0 )
NETGEN_2D_Parameters_Wing.SetOptimize( 1 )
NETGEN_2D_Parameters_Wing.SetFineness( 4 )
NETGEN_2D_Parameters_Wing.SetChordalError( -1 )
NETGEN_2D_Parameters_Wing.SetChordalErrorEnabled( 0 )
NETGEN_2D_Parameters_Wing.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_Wing.SetFuseEdges( 1 )
NETGEN_2D_Parameters_Wing.SetWorstElemMeasure( 0 )
NETGEN_2D_Parameters_Wing.SetQuadAllowed( 0 )
NETGEN_2D_Parameters_Wing.SetUseDelauney( 224 )
NETGEN_2D_Parameters_Wing.SetCheckChartBoundary( 3 )
NETGEN_2D_Parameters_Wake = smesh.CreateHypothesis('NETGEN_Parameters_2D', 'NETGENEngine')
NETGEN_2D_Parameters_Wake.SetMaxSize( 0.5 )
NETGEN_2D_Parameters_Wake.SetMinSize( 0.04 )
NETGEN_2D_Parameters_Wake.SetSecondOrder( 0 )
NETGEN_2D_Parameters_Wake.SetOptimize( 1 )
NETGEN_2D_Parameters_Wake.SetFineness( 4 )
NETGEN_2D_Parameters_Wake.SetChordalError( -1 )
NETGEN_2D_Parameters_Wake.SetChordalErrorEnabled( 0 )
NETGEN_2D_Parameters_Wake.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_Wake.SetFuseEdges( 1 )
NETGEN_2D_Parameters_Wake.SetQuadAllowed( 0 )
NETGEN_2D_Parameters_Wake.SetLocalSizeOnShape(Wake_TE, 0.001)
NETGEN_2D_Parameters_Wake.SetWorstElemMeasure( 21871 )
NETGEN_2D_Parameters_Wake.SetUseDelauney( 224 )
NETGEN_2D_Parameters_Wake.SetCheckChartBoundary( 3 )
Wake_1 = smesh.Mesh(Wake,'Wake')
status = Wake_1.AddHypothesis(NETGEN_2D_Parameters_Wake)
status = Wake_1.AddHypothesis(NETGEN_1D_2D)
isDone = Wake_1.Compute()
isDone = Fluid.Compute()
[ Inlet_1, Outlet_1, FarField_1, Walls_1, Wake_TE_2, Wing_1 ] = Fluid.GetGroups()
Domain_1 = Fluid.GroupOnGeom(Domain,'Domain',SMESH.VOLUME)
Domain_2 = Fluid.GroupOnGeom(Domain,'Domain',SMESH.NODE)
smesh.SetName(Fluid, 'Fluid')
try:
  Fluid.ExportMED( r'/media/marco-kratos-pc/datos/Examples_Compressible_Kratos/Salome_models/3-Onera/onera_m6_geometry/Mesh_Model.med', 0, 41, 1, Fluid, 1, [], '',-1, 1 )
  pass
except:
  print('ExportPartToMED() failed. Invalid file name?')
try:
  Wake_1.ExportSTL( r'/media/marco-kratos-pc/datos/Examples_Compressible_Kratos/Salome_models/3-Onera/onera_m6_geometry/Wake.stl', 1 )
  pass
except:
  print('ExportSTL() failed. Invalid file name?')
Sub_mesh_Wing = NETGEN_1D_2D.GetSubMesh()


## Set names of Mesh objects
smesh.SetName(Domain_2, 'Domain')
smesh.SetName(Sub_mesh_Wing, 'Sub-mesh_Wing')
smesh.SetName(NETGEN_1D_2D_3D.GetAlgorithm(), 'NETGEN 1D-2D-3D')
smesh.SetName(NETGEN_1D_2D.GetAlgorithm(), 'NETGEN 1D-2D')
smesh.SetName(Inlet_1, 'Inlet')
smesh.SetName(Domain_1, 'Domain')
smesh.SetName(FarField_1, 'FarField')
smesh.SetName(Outlet_1, 'Outlet')
smesh.SetName(Wing_1, 'Wing')
smesh.SetName(Walls_1, 'Walls')
smesh.SetName(Fluid.GetMesh(), 'Fluid')
smesh.SetName(Wake_1.GetMesh(), 'Wake')
smesh.SetName(Wake_TE_2, 'Wake_TE')
smesh.SetName(NETGEN_2D_Parameters_Wake, 'NETGEN 2D Parameters_Wake')
smesh.SetName(NETGEN_2D_Parameters_Wing, 'NETGEN 2D Parameters_Wing')
smesh.SetName(NETGEN_3D_Parameters_1, 'NETGEN 3D Parameters_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
