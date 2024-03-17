#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.10.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/marco/Examples_Compressible_Kratos/Salome_models/3-Onera')

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
onera_geometry = geompy.ImportIGES("/home/marco/Examples_Compressible_Kratos/Salome_models/3-Onera/onera_m6_geometry/Onera_M6_geometry.igs")
geompy.Rotate(onera_geometry, OY, 0*math.pi/180.0)
Face_1 = geompy.MakeFaceHW(40, 40, 3)
Extrusion_1 = geompy.MakePrismDXDYDZ(Face_1, 0, 20, 0)
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
[Inlet, Outlet, Wing, Walls, FarField, TE_Wing] = geompy.GetExistingSubObjects(Domain, False)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( onera_geometry, 'onera_geometry' )
geompy.addToStudy( Face_1, 'Face_1' )
geompy.addToStudy( Extrusion_1, 'Extrusion_1' )
geompy.addToStudy( Domain, 'Domain' )
geompy.addToStudyInFather( Domain, Inlet, 'Inlet' )
geompy.addToStudyInFather( Domain, Outlet, 'Outlet' )
geompy.addToStudyInFather( Domain, Wing, 'Wing' )
geompy.addToStudyInFather( Domain, Walls, 'Walls' )
geompy.addToStudyInFather( Domain, FarField, 'FarField' )
geompy.addToStudyInFather( Domain, TE_Wing, 'TE_Wing' )

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
NETGEN_3D_Parameters.SetMaxSize( 6 )
NETGEN_3D_Parameters.SetMinSize( 0.01 )
NETGEN_3D_Parameters.SetSecondOrder( 0 )
NETGEN_3D_Parameters.SetOptimize( 1 )
NETGEN_3D_Parameters.SetFineness( 5 )
NETGEN_3D_Parameters.SetGrowthRate( 0.2 )
NETGEN_3D_Parameters.SetNbSegPerEdge( 15 )
NETGEN_3D_Parameters.SetNbSegPerRadius( 2 )
# NETGEN_3D_Parameters.SetChordalError( -1 )
NETGEN_3D_Parameters.SetChordalErrorEnabled( 1 )
NETGEN_3D_Parameters.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters.SetFuseEdges( 1 )
NETGEN_3D_Parameters.SetQuadAllowed( 0 )
NETGEN_3D_Parameters.SetLocalSizeOnShape(Wing, 0.01)
Inlet_1 = Fluid.GroupOnGeom(Inlet,'Inlet',SMESH.FACE)
Outlet_1 = Fluid.GroupOnGeom(Outlet,'Outlet',SMESH.FACE)
Wing_1 = Fluid.GroupOnGeom(Wing,'Wing',SMESH.FACE)
Walls_1 = Fluid.GroupOnGeom(Walls,'Walls',SMESH.FACE)
FarField_1 = Fluid.GroupOnGeom(FarField,'FarField',SMESH.FACE)
TE_Wing_1 = Fluid.GroupOnGeom(TE_Wing,'TE_Wing',SMESH.EDGE)
isDone = Fluid.Compute()
[ Inlet_1, Outlet_1, Wing_1, Walls_1, FarField_1, TE_Wing_1 ] = Fluid.GetGroups()
Domain_1 = Fluid.GroupOnGeom(Domain,'Domain',SMESH.VOLUME)
Domain_2 = Fluid.GroupOnGeom(Domain,'Domain',SMESH.NODE)
smesh.SetName(Fluid, 'Fluid')
try:
  Fluid.ExportMED( r'/home/marco/Examples_Compressible_Kratos/Salome_models/3-Onera/Fluid.med', 0, 41, 1, Fluid, 1, [], '',-1, 1 )
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
smesh.SetName(Domain_2, 'Domain')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
