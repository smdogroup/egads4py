import numpy as np
import os
from egads4py import egads

ctx = egads.context()
ctx.setOutLevel(2)
crm = ctx.loadModel('final_surface.igs')

# Get the 'bodies' from the model
geo, oclass, mtype, lim, bodies, sense = crm.getTopology()

# Set the top and the bottom surfaces
zc = -50
x0 = [635.000, 30.057, zc]
x1 = [769.620, 30.057, zc]
x2 = [635.000, 55.033, zc]
x3 = [769.620, 55.033, zc]
x4 = [680.000, 80.010, zc]
x5 = [769.620, 80.010, zc]

X = np.array([x0, x1, x2, x3, x4, x5])
frame_edges = [[0, 1], [2, 3], [4, 5], 
               [0, 2], [1, 3], [2, 4], [3, 5]]
frame_loops = [[0, 3, 1, 4], 
               [1, 2, 3, 6]]

# Create the nodes
nodes = []
for i in range(len(X)):
    oclass = egads.NODE
    node = ctx.makeTopology(oclass, rdata=X[i])
    nodes.append(node)

# Create the lines/edges
edges = []
for e in frame_edges:
    d = X[e[1]] - X[e[0]] 
    oclass = egads.CURVE
    mtype = egads.LINE
    line = ctx.makeGeometry(oclass, mtype, rdata=[X[e[0]], d])

    topo_class = egads.EDGE
    topo_type = egads.TWONODE
    edge = ctx.makeTopology(topo_class, topo_type, geom=line,
                            children=[nodes[e[0]], nodes[e[1]]],
                            rdata=[0, np.sqrt(np.dot(d, d))])
    edges.append(edge)

# Create a series of loop
loop, nedges = ctx.makeLoop(edges)

# Create a wire body
topo_class = egads.BODY
topo_type = egads.WIREBODY
wbody = ctx.makeTopology(topo_class, topo_type, children=[loop])

# Create the extruded body
dist = 250.0
direction = [0, 0, 1]
extruded = wbody.extrude(dist, direction)

model = ctx.makeTopology(egads.MODEL, children=[extruded])
model.saveModel('extruded-model.step')

# Perform the intersection
face_list = []
for body in bodies:
    result, pairs = body.intersection(extruded)

    if result is not None:
        # Imprint the wire body onto the face
        body = body.imprintBody(pairs)

    # Get all the faces from the body
    faces = body.getBodyTopos(egads.FACE)
    for f in faces:
        face_list.append(f)

# Sew the bodies back together
print 'Sewing faces together'
print 'len(face_list) = ', len(face_list)
oml = ctx.sewFaces(face_list, toler=1e-5, manifold=False)

fname = 'test-surface-model.step'
os.system('rm %s'%(fname))

oml.saveModel(fname)
