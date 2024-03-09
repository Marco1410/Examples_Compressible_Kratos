#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.11.0 with dump python functionality
###

import os
import salome

salome.salome_init()
study = salome.myStudy
        
import salome_notebook
notebook = salome_notebook.NoteBook()

# Parameters:
AOA = 0
Domain_Length = 100.0
Domain_Width = 100.0

Airfoil_MeshSize = 5e-4
FarField_MeshSize = 10.0
Min_Airfoil_MeshSize = Airfoil_MeshSize/10.0
Deflection           = Min_Airfoil_MeshSize/10.0
salome_files_path = "salome_files"
if not os.path.exists(salome_files_path):
    os.makedirs(salome_files_path)
print (' AOA = ', AOA, 'FarField_MeshSize = ', FarField_MeshSize, 'Airfoil_MeshSize', Airfoil_MeshSize)

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS

import time as time
print(' Starting geometry ')
start_time = time.time()

geompy = geomBuilder.New()

#Create origin and axis
O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

#Create naca0012 with center in origin
Curve_UpperSurface = geompy.MakeCurveParametric("t-0.5", "0.6*(0.2969*sqrt(t) - 0.1260*t - 0.3516*t**2 + 0.2843*t**3 - 0.1036*t**4)", "0", 0, 1, 999, GEOM.Interpolation, True)
Curve_LowerSurface = geompy.MakeCurveParametric("t-0.5", "-0.6*(0.2969*sqrt(t) - 0.1260*t - 0.3516*t**2 + 0.2843*t**3 - 0.1036*t**4)", "0", 0, 1, 999, GEOM.Interpolation, True)

#Rotate around center
geompy.Rotate(Curve_UpperSurface, OZ, -AOA*math.pi/180.0)
geompy.Rotate(Curve_LowerSurface, OZ, -AOA*math.pi/180.0)

#Create face
Face_Airfoil = geompy.MakeFaceWires([Curve_UpperSurface, Curve_LowerSurface], 1)

#Create domain
Face_Domain = geompy.MakeFaceHW(Domain_Length, Domain_Width, 1)

#Cut the airfoil from the domain
Cut_Domain = geompy.MakeCutList(Face_Domain, [Face_Airfoil], True)

#Explode edges
[Edge_Inlet,Edge_WallDown,Edge_LowerSurface,Edge_UpperSurface,Edge_WallUp,Edge_Outlet] = geompy.ExtractShapes(Cut_Domain, geompy.ShapeType["EDGE"], True)
Parts_Parts_Auto1 = geompy.CreateGroup(Cut_Domain, geompy.ShapeType["FACE"])
geompy.UnionIDs(Parts_Parts_Auto1, [1])
PotentialWallCondition2D_Far_field_Auto1 = geompy.CreateGroup(Cut_Domain, geompy.ShapeType["EDGE"])
geompy.UnionIDs(PotentialWallCondition2D_Far_field_Auto1, [10, 3, 8, 6])
Body2D = geompy.CreateGroup(Cut_Domain, geompy.ShapeType["EDGE"])
geompy.UnionIDs(Body2D, [12, 15])

#Add to study
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

###
### SMESH component
###

import SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
Fluid = smesh.Mesh(Cut_Domain)

#Set NETGEN
NETGEN_1D_2D = Fluid.Triangle(algo=smeshBuilder.NETGEN_1D2D)
NETGEN_2D_Parameters_1 = NETGEN_1D_2D.Parameters()
NETGEN_2D_Parameters_1.SetMaxSize( FarField_MeshSize )
NETGEN_2D_Parameters_1.SetSecondOrder( 0 )
NETGEN_2D_Parameters_1.SetOptimize( 1 )
NETGEN_2D_Parameters_1.SetFineness( 4 )
NETGEN_2D_Parameters_1.SetMinSize( 1e-08 )
NETGEN_2D_Parameters_1.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_1.SetFuseEdges( 1 )
NETGEN_2D_Parameters_1.SetQuadAllowed( 0 )

##Set submeshes
#Body2D
Regular_1D = Fluid.Segment(geom=Body2D)
Airfoil = Regular_1D.GetSubMesh()
Adaptive_1 = Regular_1D.Adaptive(Min_Airfoil_MeshSize, Airfoil_MeshSize, Deflection)

#FarField
Regular_1D_2 = Fluid.Segment(geom=PotentialWallCondition2D_Far_field_Auto1)
FarField = Regular_1D_2.GetSubMesh()
Local_Length_1 = Regular_1D_2.LocalLength(FarField_MeshSize,None,1e-07)

Fluid_Fluid = Fluid.GroupOnGeom(Cut_Domain,'Cut_Domain',SMESH.FACE)
Fluid_Fluid.SetName( 'Fluid' )

Fluid_Airfoil = Fluid.GroupOnGeom(Body2D,'Body2D',SMESH.EDGE)
Fluid_Airfoil.SetName( 'Airfoil' )

Fluid_FarField = Fluid.GroupOnGeom(PotentialWallCondition2D_Far_field_Auto1,'PotentialWallCondition2D_Far_field_Auto1',SMESH.EDGE)
Fluid_FarField.SetName( 'FarField' )

import time as time
print(' Starting meshing ')
start_time = time.time()

#Compute mesh
isDone = Fluid.Compute()
exe_time = time.time() - start_time
print(' Mesh execution took ', str(round(exe_time, 2)), ' sec')

NumberOfNodes = Fluid.NbNodes()
NumberOfElements = Fluid.NbTriangles()
print(' Information about volume mesh:')
print(' Number of nodes       :', NumberOfNodes)
print(' Number of elements    :', NumberOfElements)

## Set names of Mesh objects
smesh.SetName(FarField, 'FarField')
smesh.SetName(NETGEN_1D_2D.GetAlgorithm(), 'NETGEN 1D-2D')
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Fluid.GetMesh(), 'Fluid')
smesh.SetName(Adaptive_1, 'Adaptive_1')
smesh.SetName(NETGEN_2D_Parameters_1, 'NETGEN 2D Parameters_1')
smesh.SetName(Local_Length_1, 'Local Length_1')
smesh.SetName(Airfoil, 'Airfoil')
smesh.SetName(Fluid_Airfoil, 'Airfoil')
smesh.SetName(Fluid_FarField, 'FarField')

Fluid.ExportMED(salome_files_path + "/model_mesh.med")

#Save salome files
salome.myStudy.SaveAs(salome_files_path + "/salome_model.hdf", study, False)

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
