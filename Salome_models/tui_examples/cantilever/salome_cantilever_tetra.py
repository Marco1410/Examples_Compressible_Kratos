#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.3.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()

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
domain = geompy.MakeBoxDXDYDZ(2, 5, 10)
[dirichlet] = geompy.SubShapes(domain, [31])
[neumann] = geompy.SubShapes(domain, [33])
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( domain, 'domain' )
geompy.addToStudyInFather( domain, dirichlet, 'dirichlet' )
geompy.addToStudyInFather( domain, neumann, 'neumann' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

domain_1 = smesh.Mesh(domain)
NETGEN_1D_2D = domain_1.Triangle(algo=smeshBuilder.NETGEN_1D2D)
NETGEN_2D_Parameters_1 = NETGEN_1D_2D.Parameters()
NETGEN_2D_Parameters_1.SetMaxSize( 1.13578 )
NETGEN_2D_Parameters_1.SetMinSize( 0.378594 )
NETGEN_2D_Parameters_1.SetSecondOrder( 0 )
NETGEN_2D_Parameters_1.SetOptimize( 1 )
NETGEN_2D_Parameters_1.SetFineness( 2 )
NETGEN_2D_Parameters_1.SetChordalError( -1 )
NETGEN_2D_Parameters_1.SetChordalErrorEnabled( 0 )
NETGEN_2D_Parameters_1.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_1.SetFuseEdges( 1 )
NETGEN_2D_Parameters_1.SetWorstElemMeasure( 0 )
NETGEN_2D_Parameters_1.SetQuadAllowed( 0 )
NETGEN_2D_Parameters_1.SetUseDelauney( 128 )
NETGEN_2D_Parameters_1.SetCheckChartBoundary( 3 )
NETGEN_3D = domain_1.Tetrahedron()
NETGEN_3D_Parameters_1 = NETGEN_3D.Parameters()
NETGEN_3D_Parameters_1.SetMaxSize( 1.13578 )
NETGEN_3D_Parameters_1.SetMinSize( 0.378594 )
NETGEN_3D_Parameters_1.SetOptimize( 1 )
NETGEN_3D_Parameters_1.SetFineness( 2 )
NETGEN_3D_Parameters_1.SetElemSizeWeight( 0 )
NETGEN_3D_Parameters_1.SetCheckOverlapping( 1 )
NETGEN_3D_Parameters_1.SetCheckChartBoundary( 3 )
isDone = domain_1.Compute()
neumann_1 = domain_1.GetSubMesh( neumann, 'Sub-mesh_2' )
group_entities_0D = domain_1.Add0DElementsToAllNodes( neumann_1, 'group_entities_0D' )
dirichlet_1 = domain_1.GetSubMesh( dirichlet, 'Sub-mesh_1' )


## Set names of Mesh objects
smesh.SetName(NETGEN_3D.GetAlgorithm(), 'NETGEN 3D')
smesh.SetName(NETGEN_1D_2D.GetAlgorithm(), 'NETGEN_1D_2D')
smesh.SetName(neumann_1, 'neumann')
smesh.SetName(dirichlet_1, 'dirichlet')
smesh.SetName(domain_1.GetMesh(), 'domain')
smesh.SetName(group_entities_0D, 'group_entities_0D')

# https://docs.salome-platform.org/latest/tui/KERNEL/kernel_salome.html
# saving the study such that it can be loaded in Salome
salome.myStudy.SaveAs("Cantilever_Tetra.hdf", False, False) # args: use_multifile, use_acsii

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()


# from here on using the plugin to create mdpa file
sys.path.append("../../../../KratosSalomePlugin/") # adding root folder of plugin to path
import create_kratos_input_tui

mesh_description_3D = { "elements"   : {"Tetra" : {"SmallDisplacementElement3D4N" : 0} } }
mesh_description_0D = { "conditions" : {"0D"    : {"PointLoadCondition3D1N"       : 0} } }

meshes = [
    create_kratos_input_tui.SalomeMesh(domain_1, mesh_description_3D, "domain"),
    create_kratos_input_tui.SalomeMesh(dirichlet_1, {}, "dirichlet"), # no elements / conditions needed
    create_kratos_input_tui.SalomeMesh(group_entities_0D, mesh_description_0D, "neumann")
]

create_kratos_input_tui.CreateMdpaFile(meshes, "cantilever_tetra")
