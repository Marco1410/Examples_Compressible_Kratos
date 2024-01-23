#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.11.0 with dump python functionality
###

import os
import salome

salome.salome_init()
study = salome.myStudy

import salome_notebook
notebook = salome_notebook.NoteBook(study)

# Parameters:
AOA = 0.0
Wing_span = 4.0
Domain_Length = 100
Domain_Height = Domain_Length
Domain_Width = 100

Airfoil_Mesh_Size = 0.002
Biggest_Airfoil_Mesh_Size = 0.002
Growth_Rate_Wing = 0.2
Growth_Rate_Domain = 0.2
Far_Field_Mesh_Size = Domain_Length/10.0
LE_Mesh_Size = Airfoil_Mesh_Size
TE_Mesh_Size = Airfoil_Mesh_Size

salome_files_path = "salome_files"
if not os.path.exists(salome_files_path):
    os.makedirs(salome_files_path)
print (' AOA = ', AOA, 'FarField_MeshSize = ', Far_Field_Mesh_Size, 'Airfoil_MeshSize', Airfoil_Mesh_Size)

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

# Create origin and axis
O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

# Create naca0012
Curve_UpperSurface_LE = geompy.MakeCurveParametric("t - 0.5", "0.0", "0.6*(0.2969*sqrt(t) - 0.1260*t - 0.3516*t**2 + 0.2843*t**3 - 0.1036*t**4)", 0, 0.5, 999, GEOM.Interpolation, True)
Curve_UpperSurface_TE = geompy.MakeCurveParametric("t - 0.5", "0.0", "0.6*(0.2969*sqrt(t) - 0.1260*t - 0.3516*t**2 + 0.2843*t**3 - 0.1036*t**4)", 0.5, 1, 999, GEOM.Interpolation, True)
Curve_LowerSurface_TE = geompy.MakeCurveParametric("t - 0.5", "0.0", "-0.6*(0.2969*sqrt(t) - 0.1260*t - 0.3516*t**2 + 0.2843*t**3 - 0.1036*t**4)", 0.5, 1, 999, GEOM.Interpolation, True)
Curve_LowerSurface_LE = geompy.MakeCurveParametric("t - 0.5", "0.0", "-0.6*(0.2969*sqrt(t) - 0.1260*t - 0.3516*t**2 + 0.2843*t**3 - 0.1036*t**4)", 0, 0.5, 999, GEOM.Interpolation, True)

# Create face
Face_Airfoil = geompy.MakeFaceWires([Curve_UpperSurface_LE, Curve_UpperSurface_TE, Curve_LowerSurface_TE, Curve_LowerSurface_LE], 1)

# Rotate around center to AOA
geompy.Rotate(Face_Airfoil, OY, AOA*math.pi/180.0)

# Extrusion of the wing
Extrusion_Wing = geompy.MakePrismVecH2Ways(Face_Airfoil, OY, Wing_span/2.0)

# Domain generation
Face_Domain = geompy.MakeFaceHW(Domain_Length, Domain_Height, 3)
Extrusion_Domain = geompy.MakePrismVecH2Ways(Face_Domain, OY, Domain_Width/2.0)

# Cut wing from the domain
Cut_Domain = geompy.MakeCutList(Extrusion_Domain, [Extrusion_Wing], True)

# Explode faces and edges
[Face_Inlet,Face_Left_Wall,Face_Left_Wing,Face_Lower_LE,Face_Upper_LE,Face_Down_Wall,Face_Top_Wall,Face_Right_Wing,Face_Lower_TE,Face_Upper_TE,Face_Right_Wall,Face_Outlet] = geompy.ExtractShapes(Cut_Domain, geompy.ShapeType["FACE"], True)

# Extruding far field
[Edge_1,Edge_2,Edge_3,Edge_4] = geompy.ExtractShapes(Face_Inlet, geompy.ShapeType["EDGE"], True)
[Obj1,Edge_6,Edge_7,Obj2] = geompy.ExtractShapes(Face_Left_Wall, geompy.ShapeType["EDGE"], True)
[Obj1,Edge_8,Edge_9,Obj2] = geompy.ExtractShapes(Face_Right_Wall, geompy.ShapeType["EDGE"], True)
[Edge_5,Edge_10,Edge_11,Edge_12] = geompy.ExtractShapes(Face_Outlet, geompy.ShapeType["EDGE"], True)

# Extruding wing
[Edge_Left_LowerLE,Edge_Left_UpperLE,Edge_Left_Lower_TE,Edge_Left_Upper_TE] = geompy.ExtractShapes(Face_Left_Wing, geompy.ShapeType["EDGE"], True)
[Edge_LE,Obj1,Obj2,Edge_Middle_Lower] = geompy.ExtractShapes(Face_Lower_LE, geompy.ShapeType["EDGE"], True)
[Obj1,Obj2,Obj3,Edge_Middle_Upper] = geompy.ExtractShapes(Face_Upper_LE, geompy.ShapeType["EDGE"], True)
[Edge_Right_LowerLE,Edge_Right_UpperLE,Edge_Right_LowerTE,Edge_Right_UpperTE] = geompy.ExtractShapes(Face_Right_Wing, geompy.ShapeType["EDGE"], True)
[Obj1,Obj2,Obj3,Edge_TE] = geompy.ExtractShapes(Face_Lower_TE, geompy.ShapeType["EDGE"], True)

# Generate wake
Vector_Wake_Direction = geompy.MakeVectorDXDYDZ(1, 0, 0)
Extrusion_Wake = geompy.MakePrismVecH(Edge_TE, Vector_Wake_Direction, Domain_Length*0.5)

# Making groups for submeshes
# Far field surface
Auto_group_for_Sub_mesh_Far_Field_Surface = geompy.CreateGroup(Cut_Domain, geompy.ShapeType["FACE"])
geompy.UnionList(Auto_group_for_Sub_mesh_Far_Field_Surface, [Face_Inlet, Face_Left_Wall, Face_Down_Wall, Face_Top_Wall, Face_Right_Wall, Face_Outlet])

# Wing surface
Auto_group_for_Sub_mesh_Wing_Surface = geompy.CreateGroup(Cut_Domain, geompy.ShapeType["FACE"])
geompy.UnionList(Auto_group_for_Sub_mesh_Wing_Surface, [Face_Left_Wing, Face_Lower_LE, Face_Upper_LE, Face_Right_Wing, Face_Lower_TE, Face_Upper_TE])

# Far field edges
Auto_group_for_Sub_mesh_Far_Field_Edges = geompy.CreateGroup(Cut_Domain, geompy.ShapeType["EDGE"])
geompy.UnionList(Auto_group_for_Sub_mesh_Far_Field_Edges, [Edge_1, Edge_2, Edge_3, Edge_4, Edge_6, Edge_7, Edge_8, Edge_9, Edge_5, Edge_10, Edge_11, Edge_12])

# TE Airfoil edges
Auto_group_for_Sub_mesh_LE_Airfoils = geompy.CreateGroup(Cut_Domain, geompy.ShapeType["EDGE"])
geompy.UnionList(Auto_group_for_Sub_mesh_LE_Airfoils, [Edge_Left_LowerLE, Edge_Left_UpperLE, Edge_Right_LowerLE, Edge_Right_UpperLE])

# LE Airfoil edges
Auto_group_for_Sub_mesh_TE_Airfoils = geompy.CreateGroup(Cut_Domain, geompy.ShapeType["EDGE"])
geompy.UnionList(Auto_group_for_Sub_mesh_TE_Airfoils, [Edge_Left_Lower_TE, Edge_Left_Upper_TE, Edge_Right_LowerTE, Edge_Right_UpperTE])

# LETE edges
Auto_group_for_Sub_mesh_LETE = geompy.CreateGroup(Cut_Domain, geompy.ShapeType["EDGE"])
geompy.UnionList(Auto_group_for_Sub_mesh_LETE, [Edge_LE, Edge_TE])

# Middle
Auto_group_for_Sub_mesh_Middle = geompy.CreateGroup(Cut_Domain, geompy.ShapeType["EDGE"])
geompy.UnionList(Auto_group_for_Sub_mesh_Middle, [Edge_Middle_Lower, Edge_Middle_Upper])

# Adding to study
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Curve_UpperSurface_LE, 'Curve_UpperSurface_LE' )
geompy.addToStudy( Curve_UpperSurface_TE, 'Curve_UpperSurface_TE' )
geompy.addToStudy( Curve_LowerSurface_TE, 'Curve_LowerSurface_TE' )
geompy.addToStudy( Curve_LowerSurface_LE, 'Curve_LowerSurface_LE' )
geompy.addToStudy( Face_Airfoil, 'Face_Airfoil' )
geompy.addToStudy( Extrusion_Wing, 'Extrusion_Wing' )
geompy.addToStudy( Face_Domain, 'Face_Domain' )
geompy.addToStudy( Extrusion_Domain, 'Extrusion_Domain' )
geompy.addToStudy( Cut_Domain, 'Cut_Domain' )
geompy.addToStudyInFather( Cut_Domain, Face_Inlet, 'Face_Inlet' )
geompy.addToStudyInFather( Cut_Domain, Face_Left_Wall, 'Face_Left_Wall' )
geompy.addToStudyInFather( Cut_Domain, Face_Left_Wing, 'Face_Left_Wing' )
geompy.addToStudyInFather( Cut_Domain, Face_Lower_LE, 'Face_Lower_LE' )
geompy.addToStudyInFather( Cut_Domain, Face_Upper_LE, 'Face_Upper_LE' )
geompy.addToStudyInFather( Cut_Domain, Face_Down_Wall, 'Face_Down_Wall' )
geompy.addToStudyInFather( Cut_Domain, Face_Top_Wall, 'Face_Top_Wall' )
geompy.addToStudyInFather( Cut_Domain, Face_Right_Wing, 'Face_Right_Wing' )
geompy.addToStudyInFather( Cut_Domain, Face_Lower_TE, 'Face_Lower_TE' )
geompy.addToStudyInFather( Cut_Domain, Face_Upper_TE, 'Face_Upper_TE' )
geompy.addToStudyInFather( Cut_Domain, Face_Right_Wall, 'Face_Right_Wall' )
geompy.addToStudyInFather( Cut_Domain, Face_Outlet, 'Face_Outlet' )
geompy.addToStudyInFather( Face_Inlet, Edge_1, 'Edge_1' )
geompy.addToStudyInFather( Face_Inlet, Edge_2, 'Edge_2' )
geompy.addToStudyInFather( Face_Inlet, Edge_3, 'Edge_3' )
geompy.addToStudyInFather( Face_Inlet, Edge_4, 'Edge_4' )
geompy.addToStudyInFather( Face_Left_Wing, Edge_Left_LowerLE, 'Edge_Left_LowerLE' )
geompy.addToStudyInFather( Face_Left_Wall, Edge_6, 'Edge_6' )
geompy.addToStudyInFather( Face_Left_Wall, Edge_7, 'Edge_7' )
geompy.addToStudyInFather( Face_Left_Wing, Edge_Left_UpperLE, 'Edge_Left_UpperLE' )
geompy.addToStudyInFather( Face_Left_Wing, Edge_Left_Lower_TE, 'Edge_Left_Lower_TE' )
geompy.addToStudyInFather( Face_Left_Wing, Edge_Left_Upper_TE, 'Edge_Left_Upper_TE' )
geompy.addToStudyInFather( Face_Lower_LE, Edge_LE, 'Edge_LE' )
geompy.addToStudyInFather( Face_Right_Wing, Edge_Right_LowerLE, 'Edge_Right_LowerLE' )
geompy.addToStudyInFather( Face_Right_Wing, Edge_Right_UpperLE, 'Edge_Right_UpperLE' )
geompy.addToStudyInFather( Face_Lower_LE, Edge_Middle_Lower, 'Edge_Middle_Lower' )
geompy.addToStudyInFather( Face_Right_Wing, Edge_Right_LowerTE, 'Edge_Right_LowerTE' )
geompy.addToStudyInFather( Face_Upper_LE, Edge_Middle_Upper, 'Edge_Middle_Upper' )
geompy.addToStudyInFather( Face_Right_Wing, Edge_Right_UpperTE, 'Edge_Right_UpperTE' )
geompy.addToStudyInFather( Face_Outlet, Edge_5, 'Edge_5' )
geompy.addToStudyInFather( Face_Right_Wall, Edge_8, 'Edge_8' )
geompy.addToStudyInFather( Face_Right_Wall, Edge_9, 'Edge_9' )
geompy.addToStudyInFather( Face_Lower_TE, Edge_TE, 'Edge_TE' )
geompy.addToStudyInFather( Face_Outlet, Edge_10, 'Edge_10' )
geompy.addToStudyInFather( Face_Outlet, Edge_11, 'Edge_11' )
geompy.addToStudyInFather( Face_Outlet, Edge_12, 'Edge_12' )
geompy.addToStudy( Extrusion_Wake, 'Extrusion_Wake' )
geompy.addToStudyInFather( Cut_Domain, Auto_group_for_Sub_mesh_Far_Field_Surface, 'Auto_group_for_Sub-mesh_Far_Field_Surface' )
geompy.addToStudyInFather( Cut_Domain, Auto_group_for_Sub_mesh_Wing_Surface, 'Auto_group_for_Sub-mesh_Wing_Surface' )
geompy.addToStudyInFather( Cut_Domain, Auto_group_for_Sub_mesh_Far_Field_Edges, 'Auto_group_for_Sub-mesh_Far_Field_Edges' )
geompy.addToStudyInFather( Cut_Domain, Auto_group_for_Sub_mesh_LE_Airfoils, 'Auto_group_for_Sub-mesh_LE_Airfoils' )
geompy.addToStudyInFather( Cut_Domain, Auto_group_for_Sub_mesh_TE_Airfoils, 'Auto_group_for_Sub-mesh_TE_Airfoils' )
geompy.addToStudyInFather( Cut_Domain, Auto_group_for_Sub_mesh_LETE, 'Auto_group_for_Sub-mesh_LETE' )
geompy.addToStudyInFather( Cut_Domain, Auto_group_for_Sub_mesh_Middle, 'Auto_group_for_Sub-mesh_Middle' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
Mesh_Domain = smesh.Mesh(Cut_Domain)

# Set NETGEN 3D
NETGEN_3D = Mesh_Domain.Tetrahedron()
NETGEN_3D_Parameters_1 = NETGEN_3D.Parameters()
NETGEN_3D_Parameters_1.SetMaxSize( Far_Field_Mesh_Size )
NETGEN_3D_Parameters_1.SetOptimize( 1 )
NETGEN_3D_Parameters_1.SetFineness( 5 )
NETGEN_3D_Parameters_1.SetGrowthRate( Growth_Rate_Domain )
NETGEN_3D_Parameters_1.SetNbSegPerEdge( 3 )
NETGEN_3D_Parameters_1.SetNbSegPerRadius( 5 )
NETGEN_3D_Parameters_1.SetMinSize( Airfoil_Mesh_Size )
NETGEN_3D_Parameters_1.SetUseSurfaceCurvature( 0 )
NETGEN_3D_Parameters_1.SetSecondOrder( 106 )
NETGEN_3D_Parameters_1.SetFuseEdges( 80 )
NETGEN_3D_Parameters_1.SetQuadAllowed( 127 )

##Set submeshes
#Far field surface
NETGEN_2D = Mesh_Domain.Triangle(algo=smeshBuilder.NETGEN_2D,geom=Auto_group_for_Sub_mesh_Far_Field_Surface)
Sub_mesh_Far_Field_Surface = NETGEN_2D.GetSubMesh()
NETGEN_2D_Parameters_FarField = NETGEN_2D.Parameters()
NETGEN_2D_Parameters_FarField.SetMaxSize( Far_Field_Mesh_Size )
NETGEN_2D_Parameters_FarField.SetOptimize( 1 )
NETGEN_2D_Parameters_FarField.SetFineness( 0 )
NETGEN_2D_Parameters_FarField.SetMinSize( Far_Field_Mesh_Size )
NETGEN_2D_Parameters_FarField.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_FarField.SetQuadAllowed( 0 )
NETGEN_2D_Parameters_FarField.SetSecondOrder( 106 )
NETGEN_2D_Parameters_FarField.SetFuseEdges( 80 )

#Wing surface
NETGEN_2D_1 = Mesh_Domain.Triangle(algo=smeshBuilder.NETGEN_2D,geom=Auto_group_for_Sub_mesh_Wing_Surface)
Sub_mesh_Wing_Surface = NETGEN_2D_1.GetSubMesh()
NETGEN_2D_Parameters_Wing = NETGEN_2D_1.Parameters()
NETGEN_2D_Parameters_Wing.SetMaxSize( Biggest_Airfoil_Mesh_Size )
NETGEN_2D_Parameters_Wing.SetOptimize( 1 )
NETGEN_2D_Parameters_Wing.SetFineness( 5 )
NETGEN_2D_Parameters_Wing.SetGrowthRate( Growth_Rate_Wing )
NETGEN_2D_Parameters_Wing.SetNbSegPerEdge( 6.92154e-310 )
NETGEN_2D_Parameters_Wing.SetNbSegPerRadius( 5.32336e-317 )
NETGEN_2D_Parameters_Wing.SetMinSize( Airfoil_Mesh_Size )
NETGEN_2D_Parameters_Wing.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_Wing.SetQuadAllowed( 0 )
NETGEN_2D_Parameters_Wing.SetSecondOrder( 106 )
NETGEN_2D_Parameters_Wing.SetFuseEdges( 80 )

#Far field edges
Regular_1D = Mesh_Domain.Segment(geom=Auto_group_for_Sub_mesh_Far_Field_Edges)
Local_Length_Far_Field = Regular_1D.LocalLength(Far_Field_Mesh_Size,None,1e-07)

#LE Airfoils
Regular_1D_1 = Mesh_Domain.Segment(geom=Auto_group_for_Sub_mesh_LE_Airfoils)
Start_and_End_Length_LE = Regular_1D_1.StartEndLength(Airfoil_Mesh_Size,Biggest_Airfoil_Mesh_Size,[])
Start_and_End_Length_LE.SetObjectEntry( 'Cut_Domain' )

#TE Airfoils
Regular_1D_2 = Mesh_Domain.Segment(geom=Auto_group_for_Sub_mesh_TE_Airfoils)
Start_and_End_Length_TE = Regular_1D_2.StartEndLength(Biggest_Airfoil_Mesh_Size,Airfoil_Mesh_Size,[])
Start_and_End_Length_TE.SetObjectEntry( 'Cut_Domain' )

#TE
Regular_1D_3 = Mesh_Domain.Segment(geom=Edge_TE)
Sub_mesh_TE = Regular_1D_3.GetSubMesh()
Local_Length_TE = Regular_1D_3.LocalLength(TE_Mesh_Size,None,1e-07)

#LE
Regular_1D_4 = Mesh_Domain.Segment(geom=Edge_LE)
Sub_mesh_LE = Regular_1D_4.GetSubMesh()
Local_Length_LE = Regular_1D_4.LocalLength(LE_Mesh_Size,None,1e-07)

#Middle
Regular_1D_4 = Mesh_Domain.Segment(geom=Auto_group_for_Sub_mesh_Middle)
Local_Length_Middle = Regular_1D_4.LocalLength(Biggest_Airfoil_Mesh_Size,None,1e-07)

Fluid_Fluid = Mesh_Domain.GroupOnGeom(Cut_Domain,'Cut_Domain',SMESH.VOLUME)
Fluid_Fluid.SetName( 'Fluid' )

Fluid_Airfoil = Mesh_Domain.GroupOnGeom(Auto_group_for_Sub_mesh_Wing_Surface,'Body2D',SMESH.FACE)
Fluid_Airfoil.SetName( 'Airfoil' )

Fluid_Airfoil = Mesh_Domain.GroupOnGeom(Edge_TE,'Body2D',SMESH.EDGE)
Fluid_Airfoil.SetName( 'Wake' )

Fluid_FarField = Mesh_Domain.GroupOnGeom(Auto_group_for_Sub_mesh_Far_Field_Surface,'PotentialWallCondition2D_Far_field_Auto1',SMESH.FACE)
Fluid_FarField.SetName( 'FarField' )

import time as time
print(' Starting meshing ')
start_time = time.time()

#Compute mesh
isDone = Mesh_Domain.Compute()
exe_time = time.time() - start_time
print(' Mesh execution took ', str(round(exe_time, 2)), ' sec')

NumberOfNodes = Mesh_Domain.NbNodes()
NumberOfElements = Mesh_Domain.NbTetras()
print(' Information about volume mesh:')
print(' Number of nodes       :', NumberOfNodes)
print(' Number of elements    :', NumberOfElements)

Sub_mesh_Far_Field_Edges = Regular_1D.GetSubMesh()
Sub_mesh_LE_Airfoils = Regular_1D_1.GetSubMesh()
Sub_mesh_TE_Airfoils = Regular_1D_2.GetSubMesh()
Sub_mesh_Middle = Regular_1D_4.GetSubMesh()

# Mesh wake and export STL
Mesh_Wake_Surface = smesh.Mesh(Extrusion_Wake)
status = Mesh_Wake_Surface.AddHypothesis(NETGEN_2D_Parameters_FarField)
NETGEN_1D_2D_2 = Mesh_Wake_Surface.Triangle(algo=smeshBuilder.NETGEN_1D2D)
isDone = Mesh_Wake_Surface.Compute()
try:
  Mesh_Wake_Surface.ExportSTL( salome_files_path + '/wake.stl', 1 )
  pass
except:
  print('ExportSTL() failed. Invalid file name?')

## Set names of Mesh objects
smesh.SetName(NETGEN_3D.GetAlgorithm(), 'NETGEN 3D')
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(NETGEN_2D.GetAlgorithm(), 'NETGEN 2D')
smesh.SetName(NETGEN_2D_Parameters_FarField, 'NETGEN 2D Parameters_FarField')
smesh.SetName(NETGEN_2D_Parameters_Wing, 'NETGEN 2D Parameters_Wing')
smesh.SetName(NETGEN_3D_Parameters_1, 'NETGEN 3D Parameters_1')
smesh.SetName(Start_and_End_Length_TE, 'Start and End Length_TE')
smesh.SetName(Local_Length_TE, 'Local Length_TE')
smesh.SetName(Local_Length_Far_Field, 'Local Length_Far_Field')
smesh.SetName(Start_and_End_Length_LE, 'Start and End Length_LE')
smesh.SetName(Local_Length_Middle, 'Local Length_Middle')
smesh.SetName(Mesh_Domain.GetMesh(), 'Mesh_Domain')
smesh.SetName(Sub_mesh_Far_Field_Edges, 'Sub-mesh_Far_Field_Edges')
smesh.SetName(Sub_mesh_Wing_Surface, 'Sub-mesh_Wing_Surface')
smesh.SetName(Sub_mesh_Far_Field_Surface, 'Sub-mesh_Far_Field_Surface')
smesh.SetName(Sub_mesh_Middle, 'Sub-mesh_Middle')
smesh.SetName(Sub_mesh_TE, 'Sub-mesh_TE')
smesh.SetName(Sub_mesh_TE_Airfoils, 'Sub-mesh_TE_Airfoils')
smesh.SetName(Sub_mesh_LE_Airfoils, 'Sub-mesh_LE_Airfoils')
smesh.SetName(Mesh_Wake_Surface, 'Mesh_Wake_Surface')

Mesh_Domain.ExportMED(salome_files_path + "/model_mesh.med")

#Save salome files
salome.myStudy.SaveAs(salome_files_path + "/salome_model.hdf", study, False)

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
