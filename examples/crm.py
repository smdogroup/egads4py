import numpy as np

from egads4py import egads

ctx = egads.context()
ctx.setOutLevel(2)
crm = ctx.loadModel('final_surface.igs')

# Get the 'bodies' from the model
geo, oclass, mtype, lim, bodies, sense = crm.getTopology()

# Create the nodes
n1 = np.array([700, 0, -25], dtype=np.float)
n2 = np.array([1200, 700, -25], dtype=np.float)
n3 = np.array([1200, 700, 200], dtype=np.float)
n4 = np.array([700, 0, 200], dtype=np.float)

oclass = egads.NODE
node1 = ctx.makeTopology(oclass, rdata=n1)
node2 = ctx.makeTopology(oclass, rdata=n2)
node3 = ctx.makeTopology(oclass, rdata=n3)
node4 = ctx.makeTopology(oclass, rdata=n4)

# Make the lines between the edges
oclass = egads.CURVE
mtype = egads.LINE
topo_class = egads.EDGE
topo_type = egads.TWONODE

d1 = n2 - n1 
l1 = ctx.makeGeometry(oclass, mtype, rdata=[n1, d1])
e1 = ctx.makeTopology(topo_class, topo_type, geom=l1,
                      children=[node1, node2],
                      rdata=[0, np.sqrt(np.dot(d1, d1))])

d2 = n3 - n2 
l2 = ctx.makeGeometry(oclass, mtype, rdata=[n2, d2])
e2 = ctx.makeTopology(topo_class, topo_type, geom=l2,
                      children=[node2, node3],
                      rdata=[0, np.sqrt(np.dot(d2, d2))])

d3 = n4 - n3
l3 = ctx.makeGeometry(oclass, mtype, rdata=[n3, d3])
e3 = ctx.makeTopology(topo_class, topo_type, geom=l3,
                      children=[node3, node4],
                      rdata=[0, np.sqrt(np.dot(d3, d3))])

d4 = n1 - n4
l4 = ctx.makeGeometry(oclass, mtype, rdata=[n4, d4])
e4 = ctx.makeTopology(topo_class, topo_type, geom=l4,
                      children=[node4, node1],
                      rdata=[0, np.sqrt(np.dot(d4, d4))])

edges = [e1, e2, e3, e4]
loop = ctx.makeLoop(edges)
mtype = egads.SREVERSE
if loop.getArea() > 0.0:
    mtype = egads.SFORWARD
face = ctx.makeFace(loop, mtype=mtype)

# Perform the intersection
face_list = []
for body in bodies:
    result, pairs = body.intersection(face)

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

oml.saveModel(fname)
