import numpy as np
import os
from egads4py import egads

def create_edges(ctx, X, frame_edges):
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

    return edges

# Set the coordinates that will be used
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

dist = 250.0
direction = [0, 0, 1]

# Set up the context
ctx = egads.context()
ctx.setOutLevel(2)
crm = ctx.loadModel('final_surface.igs')

# Get the FACEBODES from the iges file model and prepare to sew them together
geo, oclass, mtype, lim, body_surfs, sense = crm.getTopology()

# Extract the faces and sew them
all_faces = []
for body in [body_surfs[0], body_surfs[2], body_surfs[5]]:
    faces = body.getBodyTopos(egads.FACE)
    for face in faces:
        face.attributeAdd('name', egads.ATTRSTRING, 'top') 
        all_faces.append(face)    

for body in [body_surfs[1], body_surfs[3], body_surfs[6]]:
    faces = body.getBodyTopos(egads.FACE)
    for face in faces:
        face.attributeAdd('name', egads.ATTRSTRING, 'bottom') 
        all_faces.append(face)    

body_model = ctx.sewFaces(all_faces, manifold=False)

# Now 'bodies' should be a single SHEETBODY
bodies = body_model.getChildren()

# Create the edges that form the ribs/spars
edges = create_edges(ctx, X, frame_edges)
loop, nedges = ctx.makeLoop(edges)

# Create an extrude the wirebody to create the desired ribs/spars
topo_class = egads.BODY
topo_type = egads.WIREBODY
wbody = ctx.makeTopology(topo_class, topo_type, children=[loop])
extruded = wbody.extrude(dist, direction)

for face in extruded.getBodyTopos(egads.FACE):
    face.attributeAdd('name', egads.ATTRSTRING, 'box') 

new_bodies = []
for body in bodies:
    wirebody, pairs = body.intersection(extruded)
    if wirebody is not None:
        new_body = body.imprintBody(pairs)
        new_bodies.append(new_body)
        
    wirebody, pairs = extruded.intersection(body)
    if wirebody is not None:
        extruded = extruded.imprintBody(pairs)

# Now determine which new surfaces to keep and which to discard
new_bodies.append(extruded)
all_faces = []
for body in new_bodies:
    all_faces.extend(body.getBodyTopos(egads.FACE))

# Sew all the faces together
everything = ctx.sewFaces(all_faces, toler=1e-3, manifold=False)
everything.saveModel('everything.step', overwrite=True)

# Now loop over the edges and decide what should actual
wingbox_faces = []
for body in everything.getChildren():
    # Index everything based on the edge
    edge_attrs = {}
    for face in body.getBodyTopos(egads.FACE):
        attr = face.attributeRet('name')
        for edge in body.getBodyTopos(egads.EDGE, ref=face):
            index = body.indexBodyTopo(edge)
            if index in edge_attrs:
                edge_attrs[index].append(attr)
            else:
                edge_attrs[index] = [attr]

    print edge_attrs

    for face in body.getBodyTopos(egads.FACE):
        attr = face.attributeRet('name')
        if attr == 'box':
            print 'checking box edges'
            # Check if the edges have
            has_top = False
            has_bottom = False                    
            for edge in body.getBodyTopos(egads.EDGE, ref=face):
                index = body.indexBodyTopo(edge)
                for attr in edge_attrs[index]:
                    if attr == 'top':
                        has_top = True
                    elif attr == 'bottom':
                        has_bottom = True

            if has_top and has_bottom:
                print 'added box'
                wingbox_faces.append(face)
            else:
                print 'did not add box'
        elif attr == 'top' or attr == 'bottom':
            count = 0
            for edge in body.getBodyTopos(egads.EDGE, ref=face):
                index = body.indexBodyTopo(edge)
                if 'box' in edge_attrs[index]:
                    count += 1

            if count >= 2:
                wingbox_faces.append(face)
        
# Sew the faces together
oml = ctx.sewFaces(wingbox_faces, toler=1e-3, manifold=False)

fname = 'test_surface_model.step'
oml.saveModel(fname, overwrite=True)
