#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.11.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/media/marco-kratos-pc/datos/Examples_Compressible_Kratos/Salome_models/1-naca0012_2d')

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
Curve_UpperSurface = geompy.MakeCurveParametric("t-0.5", "0.6*(0.2969*sqrt(t) - 0.1260*t - 0.3516*t**2 + 0.2843*t**3 - 0.1036*t**4)", "0", 0, 1, 999, GEOM.Interpolation, True)
Curve_LowerSurface = geompy.MakeCurveParametric("t-0.5", "-0.6*(0.2969*sqrt(t) - 0.1260*t - 0.3516*t**2 + 0.2843*t**3 - 0.1036*t**4)", "0", 0, 1, 999, GEOM.Interpolation, True)
geompy.Rotate(Curve_UpperSurface, OZ, 0*math.pi/180.0)
geompy.Rotate(Curve_LowerSurface, OZ, 0*math.pi/180.0)
Face_Airfoil = geompy.MakeFaceWires([Curve_UpperSurface, Curve_LowerSurface], 1)
Face_Domain = geompy.MakeFaceHW(100, 100, 1)
Cut_Domain = geompy.MakeCutList(Face_Domain, [Face_Airfoil], True)
[Edge_Inlet,Edge_WallDown,Edge_LowerSurface,Edge_UpperSurface,Edge_WallUp,Edge_Outlet] = geompy.ExtractShapes(Cut_Domain, geompy.ShapeType["EDGE"], True)
Parts_Parts_Auto1 = geompy.CreateGroup(Cut_Domain, geompy.ShapeType["FACE"])
geompy.UnionIDs(Parts_Parts_Auto1, [1])
PotentialWallCondition2D_Far_field_Auto1 = geompy.CreateGroup(Cut_Domain, geompy.ShapeType["EDGE"])
geompy.UnionIDs(PotentialWallCondition2D_Far_field_Auto1, [10, 3, 8, 6])
Body2D = geompy.CreateGroup(Cut_Domain, geompy.ShapeType["EDGE"])
geompy.UnionIDs(Body2D, [12, 15])
Edge_Walls = geompy.CreateGroup(Cut_Domain, geompy.ShapeType["EDGE"])
geompy.UnionIDs(Edge_Walls, [6, 8])
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Curve_UpperSurface, 'Curve_UpperSurface' )
geompy.addToStudy( Curve_LowerSurface, 'Curve_LowerSurface' )
geompy.addToStudy( Face_Airfoil, 'Face_Airfoil' )
geompy.addToStudy( Face_Domain, 'Face_Domain' )
geompy.addToStudy( Cut_Domain, 'Cut_Domain' )
geompy.addToStudyInFather( Cut_Domain, Edge_Inlet, 'Edge_Inlet' )
geompy.addToStudyInFather( Cut_Domain, Edge_WallDown, 'Edge_WallDown' )
geompy.addToStudyInFather( Cut_Domain, Edge_LowerSurface, 'Edge_LowerSurface' )
geompy.addToStudyInFather( Cut_Domain, Edge_UpperSurface, 'Edge_UpperSurface' )
geompy.addToStudyInFather( Cut_Domain, Edge_WallUp, 'Edge_WallUp' )
geompy.addToStudyInFather( Cut_Domain, Edge_Outlet, 'Edge_Outlet' )
geompy.addToStudyInFather( Cut_Domain, Parts_Parts_Auto1, 'Parts_Parts_Auto1' )
geompy.addToStudyInFather( Cut_Domain, PotentialWallCondition2D_Far_field_Auto1, 'PotentialWallCondition2D_Far_field_Auto1' )
geompy.addToStudyInFather( Cut_Domain, Body2D, 'Body2D' )
geompy.addToStudyInFather( Cut_Domain, Edge_Walls, 'Edge_Walls' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Fluid = smesh.Mesh(Cut_Domain,'Fluid')
NETGEN_1D_2D = Fluid.Triangle(algo=smeshBuilder.NETGEN_1D2D)
NETGEN_2D_Parameters_1 = NETGEN_1D_2D.Parameters()
NETGEN_2D_Parameters_1.SetMaxSize( 15 )
NETGEN_2D_Parameters_1.SetSecondOrder( 0 )
NETGEN_2D_Parameters_1.SetOptimize( 1 )
NETGEN_2D_Parameters_1.SetFineness( 4 )
NETGEN_2D_Parameters_1.SetMinSize( 1e-08 )
NETGEN_2D_Parameters_1.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_1.SetFuseEdges( 1 )
NETGEN_2D_Parameters_1.SetQuadAllowed( 0 )
Regular_1D = Fluid.Segment(geom=Body2D)
Adaptive_1 = Regular_1D.Adaptive(5e-05,0.0005,5e-06)
Regular_1D_1 = Fluid.Segment(geom=PotentialWallCondition2D_Far_field_Auto1)
Local_Length_1 = Regular_1D_1.LocalLength(15,None,1e-07)
Fluid_1 = Fluid.GroupOnGeom(Cut_Domain,'Cut_Domain',SMESH.FACE)
Fluid_1.SetName( 'Fluid' )
Airfoil = Fluid.GroupOnGeom(Body2D,'Body2D',SMESH.EDGE)
Airfoil.SetName( 'Airfoil' )
FarField = Fluid.GroupOnGeom(PotentialWallCondition2D_Far_field_Auto1,'PotentialWallCondition2D_Far_field_Auto1',SMESH.EDGE)
FarField.SetName( 'FarField' )
isDone = Fluid.Compute()
smesh.SetName(Fluid, 'Fluid')
try:
  Fluid.ExportMED(r'salome_files/model_mesh.med',auto_groups=0,version=-1,overwrite=1,meshPart=None,autoDimension=1)
  pass
except:
  print('ExportMED() failed. Invalid file name?')
Inlet = Fluid.GroupOnGeom(Edge_Inlet,'Edge_Inlet',SMESH.EDGE)
Outlet = Fluid.GroupOnGeom(Edge_Outlet,'Edge_Outlet',SMESH.EDGE)
Inlet.SetName( 'Inlet' )
Outlet.SetName( 'Outlet' )
Edge_Walls_1 = Fluid.GroupOnGeom(Edge_Walls,'Edge_Walls',SMESH.EDGE)
smesh.SetName(Fluid, 'Fluid')
try:
  Fluid.ExportMED( r'/media/marco-kratos-pc/datos/Examples_Compressible_Kratos/Salome_models/1-naca0012_2d/salome_files/model_mesh.med', 0, 41, 1, Fluid, 1, [], '',-1, 1 )
  pass
except:
  print('ExportPartToMED() failed. Invalid file name?')
[ Fluid_1, Airfoil, FarField, Inlet, Outlet, Edge_Walls_1 ] = Fluid.GetGroups()
Airfoil_1 = Regular_1D.GetSubMesh()
FarField_1 = Regular_1D_1.GetSubMesh()


## Set names of Mesh objects
smesh.SetName(NETGEN_1D_2D.GetAlgorithm(), 'NETGEN 1D-2D')
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Adaptive_1, 'Adaptive_1')
smesh.SetName(Local_Length_1, 'Local Length_1')
smesh.SetName(NETGEN_2D_Parameters_1, 'NETGEN 2D Parameters_1')
smesh.SetName(Fluid_1, 'Fluid')
smesh.SetName(Fluid.GetMesh(), 'Fluid')
smesh.SetName(FarField_1, 'FarField')
smesh.SetName(Airfoil_1, 'Airfoil')
smesh.SetName(Airfoil, 'Airfoil')
smesh.SetName(Inlet, 'Inlet')
smesh.SetName(FarField, 'FarField')
smesh.SetName(Edge_Walls_1, 'Edge_Walls')
smesh.SetName(Outlet, 'Outlet')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
