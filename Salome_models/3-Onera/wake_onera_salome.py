#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.10.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()

Wing_angle              = 0.0
Wake_length             = 5.0
Wake_angle              = 3.06
Wake_max_mesh_size      = 1.0
Wake_min_mesh_size      = 0.001
TE_Wing_mesh_size       = 0.001

Dir                     = "/media/marco-kratos-pc/datos/Examples_Compressible_Kratos/Salome_models/3-Onera"
onera_geometry_igs_path = Dir + "/onera_m6_geometry/Onera_M6_geometry.igs"
Wake_output_name        = Dir + "/Wake.stl"

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
onera_geometry = geompy.ImportIGES(onera_geometry_igs_path)
geompy.Rotate(onera_geometry, OY, Wing_angle*math.pi/180.0)
TE_wing_base = geompy.CreateGroup(onera_geometry, geompy.ShapeType["EDGE"])
geompy.UnionIDs(TE_wing_base, [16])
wake_stl = geompy.MakePrismDXDYDZ(TE_wing_base, Wake_length, 0, 0)
geompy.Rotate(wake_stl, TE_wing_base, Wake_angle*math.pi/180.0)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( onera_geometry, 'onera_geometry' )
geompy.addToStudyInFather( onera_geometry, TE_wing_base, 'TE_wing_base' )
geompy.addToStudy( wake_stl, 'wake_stl' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Wake = smesh.Mesh(wake_stl,'Wake')
NETGEN_1D_2D = Wake.Triangle(algo=smeshBuilder.NETGEN_1D2D)
NETGEN_2D_Parameters = NETGEN_1D_2D.Parameters()
NETGEN_2D_Parameters.SetMaxSize( Wake_max_mesh_size )
NETGEN_2D_Parameters.SetMinSize( Wake_min_mesh_size )
NETGEN_2D_Parameters.SetSecondOrder( 0 )
NETGEN_2D_Parameters.SetOptimize( 1 )
NETGEN_2D_Parameters.SetFineness( 2 )
NETGEN_2D_Parameters.SetChordalError( -1 )
NETGEN_2D_Parameters.SetChordalErrorEnabled( 0 )
NETGEN_2D_Parameters.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters.SetFuseEdges( 1 )
NETGEN_2D_Parameters.SetWorstElemMeasure( 0 )
NETGEN_2D_Parameters.SetQuadAllowed( 0 )
NETGEN_2D_Parameters.SetLocalSizeOnShape(TE_wing_base, TE_Wing_mesh_size)
NETGEN_2D_Parameters.SetUseDelauney( 128 )
NETGEN_2D_Parameters.SetCheckChartBoundary( 3 )
isDone = Wake.Compute()
try:
  Wake.ExportSTL( Wake_output_name, 1 )
  pass
except:
  print('ExportSTL() failed. Invalid file name?')


## Set names of Mesh objects
smesh.SetName(NETGEN_1D_2D.GetAlgorithm(), 'NETGEN 1D-2D')
smesh.SetName(NETGEN_2D_Parameters, 'NETGEN 2D Parameters')
smesh.SetName(Wake.GetMesh(), 'Wake')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
