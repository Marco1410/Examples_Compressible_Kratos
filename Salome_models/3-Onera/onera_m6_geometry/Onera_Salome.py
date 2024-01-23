#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.11.0 with dump python functionality
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

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
