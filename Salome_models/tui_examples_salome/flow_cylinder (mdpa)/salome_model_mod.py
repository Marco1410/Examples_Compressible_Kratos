#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.11.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/media/marco-kratos-pc/datos/Examples_Compressible_Kratos/Salome_models/tui_examples_salome/flow_cylinder')

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
Vertex_1 = geompy.MakeVertex(0, 0, 0)
Vertex_2 = geompy.MakeVertex(7, 0, 0)
Vertex_3 = geompy.MakeVertex(7, 2, 0)
Vertex_4 = geompy.MakeVertex(0, 2, 0)
Vertex_5 = geompy.MakeVertex(2.5, 1, 0)
Vertex_6 = geompy.MakeVertex(1.5, 0.5, 0)
Vertex_7 = geompy.MakeVertex(1.5, 1.5, 0)
Vertex_8 = geompy.MakeVertex(5, 0.5, 0)
Vertex_9 = geompy.MakeVertex(5, 1.5, 0)
Circle_1 = geompy.MakeCircle(Vertex_5, None, 0.15)
outer_boundary = geompy.MakePolyline([Vertex_1, Vertex_4, Vertex_3, Vertex_2], True)
inner_boundary = geompy.MakePolyline([Vertex_6, Vertex_7, Vertex_9, Vertex_8], True)
Face_1 = geompy.MakeFaceWires([Circle_1, outer_boundary], 1)
domain = geompy.MakePartition([Face_1], [inner_boundary], [], [], geompy.ShapeType["FACE"], 0, [], 0)
[Plane_2] = geompy.SubShapes(domain, [21])
[Edge_4,Edge_1,cyl_boundary] = geompy.SubShapes(domain, [11, 4, 24])
walls = geompy.CreateGroup(domain, geompy.ShapeType["EDGE"])
geompy.UnionIDs(walls, [7, 9])
boundary_inner = geompy.CreateGroup(domain, geompy.ShapeType["EDGE"])
geompy.UnionIDs(boundary_inner, [20, 18, 13, 16])
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Vertex_1, 'Vertex_1' )
geompy.addToStudy( Vertex_2, 'Vertex_2' )
geompy.addToStudy( Vertex_3, 'Vertex_3' )
geompy.addToStudy( Vertex_4, 'Vertex_4' )
geompy.addToStudy( Vertex_5, 'Vertex_5' )
geompy.addToStudy( Vertex_6, 'Vertex_6' )
geompy.addToStudy( Vertex_7, 'Vertex_7' )
geompy.addToStudy( Vertex_8, 'Vertex_8' )
geompy.addToStudy( Vertex_9, 'Vertex_9' )
geompy.addToStudy( Circle_1, 'Circle_1' )
geompy.addToStudy( outer_boundary, 'outer_boundary' )
geompy.addToStudy( inner_boundary, 'inner_boundary' )
geompy.addToStudy( Face_1, 'Face_1' )
geompy.addToStudy( domain, 'domain' )
geompy.addToStudyInFather( domain, Plane_2, 'Plane_2' )
geompy.addToStudyInFather( domain, Edge_4, 'Edge_4' )
geompy.addToStudyInFather( domain, Edge_1, 'Edge_1' )
geompy.addToStudyInFather( domain, cyl_boundary, 'cyl_boundary' )
geompy.addToStudyInFather( domain, walls, 'walls' )
geompy.addToStudyInFather( domain, boundary_inner, 'boundary_inner' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

domain_1 = smesh.Mesh(domain,'domain')
Regular_1D = domain_1.Segment()
Max_Size_domain = Regular_1D.MaxSize(2)
NETGEN_2D = domain_1.Triangle(algo=smeshBuilder.NETGEN_2D)
Regular_1D_1 = domain_1.Segment(geom=Plane_2)
Max_Size_domain_inner = Regular_1D_1.MaxSize(0.2)
NETGEN_2D_1 = domain_1.Triangle(algo=smeshBuilder.NETGEN_2D,geom=Plane_2)
Regular_1D_2 = domain_1.Segment(geom=Edge_1)
Local_Length_outter_boundary = Regular_1D_2.LocalLength(0.1,None,1e-07)
Regular_1D_3 = domain_1.Segment(geom=Edge_4)
status = domain_1.AddHypothesis(Local_Length_outter_boundary,Edge_4)
Regular_1D_4 = domain_1.Segment(geom=walls)
status = domain_1.AddHypothesis(Local_Length_outter_boundary,walls)
Regular_1D_5 = domain_1.Segment(geom=cyl_boundary)
Local_Length_cyl = Regular_1D_5.LocalLength(0.1,None,1e-07)
Regular_1D_6 = domain_1.Segment(geom=boundary_inner)
NETGEN_2D_2 = domain_1.Triangle(algo=smeshBuilder.NETGEN_2D)
isDone = domain_1.Compute()
domain_2 = domain_1.GroupOnGeom(domain,'domain',SMESH.FACE)
domain_3 = domain_1.GroupOnGeom(domain,'domain',SMESH.NODE)
smesh.SetName(domain_1, 'domain')
try:
  domain_1.ExportMED( r'/media/marco-kratos-pc/datos/Examples_Compressible_Kratos/Salome_models/tui_examples_salome/flow_cylinder/domain.med', 1, 41, 1, domain_1, 1, [], 'efv',-1, 1 )
  pass
except:
  print('ExportPartToMED() failed. Invalid file name?')
mesh_domain_inner = Regular_1D_1.GetSubMesh()
mesh_inlet = Regular_1D_2.GetSubMesh()
mesh_outlet = Regular_1D_3.GetSubMesh()
mesh_walls = Regular_1D_4.GetSubMesh()
mesh_cyl_boundary = Regular_1D_5.GetSubMesh()
Regular_1D_7 = Regular_1D_6.GetSubMesh()


## Set names of Mesh objects
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(NETGEN_2D.GetAlgorithm(), 'NETGEN 2D')
smesh.SetName(Max_Size_domain_inner, 'Max Size_domain_inner')
smesh.SetName(Local_Length_outter_boundary, 'Local Length_outter_boundary')
smesh.SetName(Max_Size_domain, 'Max Size_domain')
smesh.SetName(domain_2, 'domain')
smesh.SetName(Local_Length_cyl, 'Local Length_cyl')
smesh.SetName(mesh_domain_inner, 'mesh_domain_inner')
smesh.SetName(domain_1.GetMesh(), 'domain')
smesh.SetName(Regular_1D_7, 'Regular_1D')
smesh.SetName(mesh_walls, 'mesh_walls')
smesh.SetName(domain_3, 'domain')
smesh.SetName(mesh_cyl_boundary, 'mesh_cyl_boundary')
smesh.SetName(mesh_outlet, 'mesh_outlet')
smesh.SetName(mesh_inlet, 'mesh_inlet')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
