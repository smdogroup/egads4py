from egads4py import egads
from mpi4py import MPI
from tmr import TMR
import numpy as np
import argparse
import os

'''
Example to generate a 2D rectangular plate  with curved edges using egads4py
'''

# Create the plate geometry
l = 4.0 # length of the plate
h = 1.5 # height of the plate

# Create the egads context
ctx = egads.context()
ctx.setOutLevel(1)

x1 = [0.0, 0.0, 0.0]
x2 = [l,   0.0, 0.0]
x3 = [l,   h,   0.0]
x4 = [0.0, h,   0.0]
d1 = np.array(x2) - np.array(x1)
d2 = np.array(x3) - np.array(x2)
d3 = np.array(x4) - np.array(x3)
d4 = np.array(x1) - np.array(x4)

# Make the nodes
n1 = ctx.makeTopology(egads.NODE, rdata=x1)
n2 = ctx.makeTopology(egads.NODE, rdata=x2)
n3 = ctx.makeTopology(egads.NODE, rdata=x3)
n4 = ctx.makeTopology(egads.NODE, rdata=x4)

# Make the curves
l1 = ctx.makeGeometry(egads.CURVE, mtype=egads.LINE,
                      rdata=[x1, d1])
l2 = ctx.makeGeometry(egads.CURVE, mtype=egads.LINE,
                      rdata=[x2, d2])
l3 = ctx.makeGeometry(egads.CURVE, mtype=egads.LINE,
                      rdata=[x3, d3])
l4 = ctx.makeGeometry(egads.CURVE, mtype=egads.LINE,
                      rdata=[x4, d4])

# Make the edges
e1 = ctx.makeTopology(egads.EDGE, mtype=egads.TWONODE, geom=l1,
                      children=[n1, n2],
                      rdata=[0, np.sqrt(np.dot(d1, d1))])
e2 = ctx.makeTopology(egads.EDGE, mtype=egads.TWONODE, geom=l2,
                      children=[n2, n3],
                      rdata=[0, np.sqrt(np.dot(d2, d2))])
e3 = ctx.makeTopology(egads.EDGE, mtype=egads.TWONODE, geom=l3,
                      children=[n3, n4],
                      rdata=[0, np.sqrt(np.dot(d3, d3))])
e4 = ctx.makeTopology(egads.EDGE, mtype=egads.TWONODE, geom=l4,
                      children=[n4, n1],
                      rdata=[0, np.sqrt(np.dot(d4, d4))])

# Make the edge loop
el1, nloop_edges = ctx.makeLoop([e1, e2, e3, e4])

# Make the face
f1 = ctx.makeFace(el1, egads.SFORWARD)

# Create the model
s1 = ctx.makeTopology(egads.SHELL, egads.OPEN, children=[f1])
panel_body = ctx.makeTopology(egads.BODY, egads.SHEETBODY,
                              children=[s1])
panel_model = ctx.makeTopology(egads.MODEL,
                               children=[panel_body])

# Save the model
panel_model.saveModel('rect_plate.step', overwrite=True)
