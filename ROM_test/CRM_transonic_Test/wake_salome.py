#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.11.0 with dump python functionality
###

import os
import numpy as np
import sys
import salome# type: ignore

salome.salome_init()
study = salome.myStudy
import salome_notebook# type: ignore
notebook = salome_notebook.NoteBook()

PlaneAngle = 0.0
########################
Wake_angles        = np.loadtxt('wake_angles.dat', usecols=(0,))
Wake_length        = 125.0
Wake_max_mesh_size = 50.0
Wake_min_mesh_size = 0.1
TE_Plane_mesh_size = 0.1

Dir                   = os.getcwd()
crm_geometry_igs_path = Dir + "/Mesh/crm_geometry/nasa_crm_original_TE_00.igs"

if not os.path.exists('SalomeFiles'):
  os.mkdir('SalomeFiles')

###
### GEOM component
###

import GEOM # type: ignore
from salome.geom import geomBuilder# type: ignore
import math
import SALOMEDS # type: ignore


geompy = geomBuilder.New()

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
SOLID_1_ = geompy.ImportIGES(crm_geometry_igs_path, True)
CRM = geompy.MakeScaleTransform(SOLID_1_, None, 0.025)
geompy.Rotate(CRM, OY, PlaneAngle*math.pi/180.0)

Plane_TE_Sup = geompy.CreateGroup(CRM, geompy.ShapeType["EDGE"])
geompy.UnionIDs(Plane_TE_Sup, [167, 181, 211])

Wing_Root_Point_Sup = geompy.CreateGroup(CRM, geompy.ShapeType["VERTEX"])
geompy.UnionIDs(Wing_Root_Point_Sup, [152])
Tail_Root_Point_Sup = geompy.CreateGroup(CRM, geompy.ShapeType["VERTEX"])
geompy.UnionIDs(Tail_Root_Point_Sup, [120])

Wing_Root_Line_Sup = geompy.MakePrismDXDYDZ(Wing_Root_Point_Sup, 0, -4, 0)
Wing_Root_Tail_Sup = geompy.MakePrismDXDYDZ(Tail_Root_Point_Sup, 0, -2, 0)

Plane_TE_Lines = geompy.MakeFuseList([Wing_Root_Line_Sup, Wing_Root_Tail_Sup, Plane_TE_Sup], True, True)

for i, Wake_angle in enumerate(Wake_angles):

  Wake_output_name = f'{Dir}/SalomeFiles/Wake_{Wake_angle}.stl'

  Wake_sup = geompy.MakePrismDXDYDZ(Plane_TE_Lines, 500, 0, 500*math.tan(Wake_angle*math.pi/180))

  geompy.addToStudy( O, 'O' )
  geompy.addToStudy( OX, 'OX' )
  geompy.addToStudy( OY, 'OY' )
  geompy.addToStudy( OZ, 'OZ' )
  geompy.addToStudy( SOLID_1_, 'SOLID(1)' )
  geompy.addToStudy( CRM, 'CRM' )
  geompy.addToStudyInFather( CRM, Plane_TE_Sup, 'Plane_TE_Sup' )
  geompy.addToStudyInFather( CRM, Wing_Root_Point_Sup, 'Wing_Root_Point_Sup' )
  geompy.addToStudyInFather( CRM, Tail_Root_Point_Sup, 'Tail_Root_Point_Sup' )
  geompy.addToStudy( Wing_Root_Line_Sup, 'Wing_Root_Line_Sup' )
  geompy.addToStudy( Wing_Root_Tail_Sup, 'Wing_Root_Tail_Sup' )
  geompy.addToStudy( Plane_TE_Lines, 'Plane_TE_Lines' )
  geompy.addToStudy( Wake_sup, 'Wake_sup' )

  ###
  ### SMESH component
  ###

  import  SMESH, SALOMEDS # type: ignore
  from salome.smesh import smeshBuilder # type: ignore

  smesh = smeshBuilder.New()
  #smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                  # multiples meshes built in parallel, complex and numerous mesh edition (performance)

  Wake = smesh.Mesh(Wake_sup,'Wake')
  NETGEN_1D_2D_1 = Wake.Triangle(algo=smeshBuilder.NETGEN_1D2D)
  NETGEN_2D_Parameters_Wake = NETGEN_1D_2D_1.Parameters()
  NETGEN_2D_Parameters_Wake.SetMaxSize( Wake_max_mesh_size )
  NETGEN_2D_Parameters_Wake.SetMinSize( Wake_min_mesh_size )
  NETGEN_2D_Parameters_Wake.SetSecondOrder( 0 )
  NETGEN_2D_Parameters_Wake.SetOptimize( 1 )
  NETGEN_2D_Parameters_Wake.SetFineness( 2 )
  NETGEN_2D_Parameters_Wake.SetChordalError( -1 )
  NETGEN_2D_Parameters_Wake.SetChordalErrorEnabled( 0 )
  NETGEN_2D_Parameters_Wake.SetUseSurfaceCurvature( 1 )
  NETGEN_2D_Parameters_Wake.SetFuseEdges( 1 )
  NETGEN_2D_Parameters_Wake.SetQuadAllowed( 0 )
  NETGEN_2D_Parameters_Wake.SetLocalSizeOnShape(Plane_TE_Lines, Wake_min_mesh_size)
  NETGEN_2D_Parameters_Wake.SetWorstElemMeasure( 24392 )
  NETGEN_2D_Parameters_Wake.SetUseDelauney( 96 )
  NETGEN_2D_Parameters_Wake.SetCheckChartBoundary( 3 )

  isDone = Wake.Compute()

  try:
    Wake.ExportSTL( Wake_output_name, 1 )
    pass
  except:
    print('ExportSTL() failed. Invalid file name?')

  ## Set names of Mesh objects
  smesh.SetName(Wake.GetMesh(), 'Wake')
  smesh.SetName(NETGEN_2D_Parameters_Wake, 'NETGEN 2D Parameters_Wake')

  #Save salome files
  salome.myStudy.SaveAs(f'{Dir}/SalomeFiles/salome_wake_model_{Wake_angle}.hdf', study, False)

  ##########################################################################

  print(' Information about wake mesh:')
  print(' Number of nodes       :', Wake.NbNodes())
  print(' Number of elements    :', Wake.NbTriangles())

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()