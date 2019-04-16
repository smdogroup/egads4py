from mpi4py import MPI
from tmr import TMR
from egads4py import egads
import numpy as np
from coincident_geom import *

ctx = egads.context()

# Create the bodies that are 
x = [0, 0, 0]
d = [0.5, 0.5, 0.5]
b1 = ctx.makeSolidBody(egads.BOX, rdata=[x, d])

x = [0, 0, 0.5]
d = [0.5, 0.5, 0.5]
b2 = ctx.makeSolidBody(egads.BOX, rdata=[x, d])

# Write out to a STEP file
m1 = ctx.makeTopology(egads.MODEL, children=[b1])
m2 = ctx.makeTopology(egads.MODEL, children=[b2])
m1.saveModel('box1.egads', overwrite=True)
m2.saveModel('box2.egads', overwrite=True)

comm = MPI.COMM_WORLD
htarget = 0.05

geo1 = TMR.LoadModel('box1.egads', print_lev=1)
geo2 = TMR.LoadModel('box2.egads', print_lev=1)

findMatchingFaces(geo1, geo2)

verts = []
edges = []
faces = []
vols = []
for geo in [geo1, geo2]:
    f = geo.getFaces()
    v = geo.getVolumes()
    f[0].setSource(v[0], f[1])

    verts.extend(geo.getVertices())
    edges.extend(geo.getEdges())
    faces.extend(f)
    vols.extend(geo.getVolumes())

geo = TMR.Model(verts, edges, faces, vols)

# Create the new mesh
mesh = TMR.Mesh(comm, geo)

# Set the meshing options
opts = TMR.MeshOptions()
opts.write_mesh_quality_histogram = 1
opts.triangularize_print_iter = 50000

# Create the surface mesh
mesh.mesh(htarget, opts)

# Write the surface mesh to a file
mesh.writeToVTK('output.vtk', 'hex')

