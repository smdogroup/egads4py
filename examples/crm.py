import numpy as np
import os
import sys
from egads4py import egads
from dcel import dcel

def create_faces(ctx, X, frame_edges, nodes=None, edges=None, faces=None):
    # Create the nodes
    nx = len(X)
    ne = len(frame_edges)

    # Create the edges if not created already
    if nodes is None:
        nodes = 2*nx*[None]
    if edges is None:
        edges = (nx + 2*ne)*[None]
    if faces is None:
        faces = ne*[None]

    # Create the bottom nodes
    for i in range(nx):
        if nodes[i] is None:
            oclass = egads.NODE
            node = ctx.makeTopology(oclass, rdata=X[i,0])
            nodes[i] = node

    # Create the top nodes
    for i in range(nx):
        if nodes[i + nx] is None:
            oclass = egads.NODE
            node = ctx.makeTopology(oclass, rdata=X[i,1])
            nodes[i + nx] = node

    # Create the bottom edges
    for index, e in enumerate(frame_edges):
        if edges[index] is None:
            d = X[e[1],0] - X[e[0],0] 
            oclass = egads.CURVE
            mtype = egads.LINE
            line = ctx.makeGeometry(oclass, mtype, rdata=[X[e[0],0], d])
            
            topo_class = egads.EDGE
            topo_type = egads.TWONODE
            edge = ctx.makeTopology(topo_class, topo_type, geom=line,
                                    children=[nodes[e[0]], nodes[e[1]]],
                                    rdata=[0, np.sqrt(np.dot(d, d))])
            edge.attributeAdd('name', egads.ATTRSTRING, 'bottom')
            edges[index] = edge

    # Create the top edges
    for index, e in enumerate(frame_edges):
        if edges[index + ne] is None:
            d = X[e[1],1] - X[e[0],1] 
            oclass = egads.CURVE
            mtype = egads.LINE
            line = ctx.makeGeometry(oclass, mtype, rdata=[X[e[0],1], d])

            topo_class = egads.EDGE
            topo_type = egads.TWONODE
            edge = ctx.makeTopology(topo_class, topo_type, geom=line,
                                    children=[nodes[e[0]+nx], 
                                              nodes[e[1]+nx]],
                                    rdata=[0, np.sqrt(np.dot(d, d))])
            edge.attributeAdd('name', egads.ATTRSTRING, 'top')
            edges[index + ne] = edge

    # Create the bottom to top edges
    for i in range(nx):
        if edges[i + 2*ne] is None:
            d = X[i,1] - X[i,0]
            oclass = egads.CURVE
            mtype = egads.LINE
            line = ctx.makeGeometry(oclass, mtype, rdata=[X[i,0], d])

            topo_class = egads.EDGE
            topo_type = egads.TWONODE
            edge = ctx.makeTopology(topo_class, topo_type, geom=line,
                                    children=[nodes[i], nodes[nx+i]],
                                    rdata=[0, np.sqrt(np.dot(d, d))])
            edges[i + 2*ne] = edge

    # Create the faces
    for i, e in enumerate(frame_edges):
        if faces[i] is None:
            elist = [edges[i], edges[2*ne+e[1]], 
                     edges[ne+i], edges[2*ne+e[0]]]
            senses = [1, 1, -1, -1]
            loop = ctx.makeTopology(egads.LOOP, egads.CLOSED, 
                                    children=elist, sens=senses)
            mtype = egads.SFORWARD
            face = ctx.makeFace(loop, mtype=mtype)
            face.attributeAdd('name', egads.ATTRSTRING, 'box%d'%(i))
            faces[i] = face

    return faces

def load_oml_model(ctx, leList, teList, omlfile):
    # Load the model from a step file
    crm = ctx.loadModel(omlfile)

    # Get the SHELLBODY or FACEBODY from the model and prepare 
    # to sew them together
    body = crm.getChildren()[0]

    # Extract the faces 
    surf_faces = []
    for face in body.getBodyTopos(egads.FACE):
        surf_faces.append(face)

    # Go through the edge list and find whether the edge is in front
    # of the leading edge or behind the trailing edge
    for edge in body.getBodyTopos(egads.EDGE):
        # Get the range of parameter values
        r, periodic = edge.getRange()

        # Evalaute the 
        types = [None, None]
        for k in range(2):
            x, xt, xtt = edge.evaluate(r[k])

            # Evaluate the 
            y = x[1]

            # Compute the leading and trailing edge locations
            xle = 0.0
            xte = 0.0
            if y < leList[0,1]:
                xle = leList[0,0]
            elif y > leList[-1,1]:
                xle = leList[-1,0]
            else:
                for j in range(len(leList)-1):
                    if (y >= leList[j,1] and y < leList[j+1,1]):
                        xi = (y - leList[j,1])/(leList[j+1,1] - leList[j,1])
                        xle = (1.0 - xi)*leList[j,0] + xi*leList[j+1,0]
                        break

            if y < teList[0,1]:
                xte = teList[0,0]
            elif y > teList[-1,1]:
                xte = teList[-1,0]
            else:
                for j in range(len(teList)-1):
                    if (y >= teList[j,1] and y < teList[j+1,1]):
                        xi = (y - teList[j,1])/(teList[j+1,1] - teList[j,1])
                        xte = (1.0 - xi)*teList[j,0] + xi*teList[j+1,0]
                        break

            if x[0] > xle:
                types[k] = 'leading'
            elif x[0] < xte:
                types[k] = 'trailing'

        if types[0] == types[1]:
            edge.attributeAdd('name', egads.ATTRSTRING, types[0])

    return body

def create_crm_model(ctx, body, ribspars, filename='ucrm.step'):
    # Keep track of the faces that should be in the final wing-box
    wingbox_faces = []

    # Set the faces
    print 'Intersecting iges surface body...'
    wirebody, body_pairs = body.intersection(ribspars)
    if wirebody is not None:
        print 'Impriting iges surface body...'
        new_body = body.imprintBody(body_pairs)

    print 'Intersecting ribspar body...'
    wirebody, ribspar_pairs = ribspars.intersection(body)
    if wirebody is not None:
        print 'Imprinting ribspar body...'
        ribspars = ribspars.imprintBody(ribspar_pairs)

    # Go through and discard surfaces not cut by bodies
    for face in new_body.getBodyTopos(egads.FACE):
        add_edge = True
        for edge in new_body.getBodyTopos(egads.EDGE, ref=face):
            if 'leading' == edge.attributeRet('name'):
                add_edge = False
            elif 'trailing' == edge.attributeRet('name'):
                add_edge = False
        if add_edge:
            wingbox_faces.append(face)

    # Go through the ribspar body and remove the top/bottom 
    # pieces of the ribs and spars 
    for face in ribspars.getBodyTopos(egads.FACE):
        has_top = False
        has_bottom = False
        for edge in ribspars.getBodyTopos(egads.EDGE, ref=face):
            # Find the edge created by the intersection and compare
            # the face name to see if it is cut by both the top
            # and bottom surfaces
            if 'top' == edge.attributeRet('name'):
                has_top = True
            elif 'bottom' == edge.attributeRet('name'):
                has_bottom = True

        if has_top is False and has_bottom is False:
            wingbox_faces.append(face)
    
    # Sew the faces together
    print 'Sewing it all together...'
    oml = ctx.sewFaces(wingbox_faces, toler=1e-4, manifold=False)
    oml.saveModel(filename, overwrite=True)

    # # Create the model
    # shell = ctx.makeTopology(egads.SHELL, egads.OPEN, children=wingbox_faces)
    # body = ctx.makeTopology(egads.BODY, egads.SHEETBODY, children=[shell])
    # oml = ctx.makeTopology(egads.MODEL, children=[body])
    # oml.saveModel(filename, overwrite=True)

    return

def compute_ribspar_edges(leList, teList, nrib1=5, nrib2=44):
    # Compute and normalize the rib direction that is norma
    # trailing edge direction
    xspar = np.array([teList[1,0], teList[1,1]])
    dspar = np.array([teList[2,0] - teList[1,0],
                      teList[2,1] - teList[1,1]])
    dspar = dspar/np.sqrt(np.dot(dspar, dspar))
    
    # Compute the direction of all ribs (normal to the te spar)
    drib = np.array([dspar[1], -dspar[0]])

    # Compute the first point and direction of the root rib
    xroot = np.array([leList[1,0], leList[1,1]])
    droot = np.array([teList[1,0] - leList[1,0],
                      teList[1,1] - leList[1,1]])

    # Create the connectivity
    conn = [[0, 1, 4, 3], 
            [1, 2, 5, 4]]
    X = [[teList[0,0], teList[0,1]],
         [teList[1,0], teList[1,1]],
         [teList[2,0], teList[2,1]],
         [leList[0,0], leList[0,1]],
         [leList[1,0], leList[1,1]],
         [leList[2,0], leList[2,1]]]

    # Create the dcel object
    d = dcel(X, conn)

    # Find the intersection with the through-body box
    x = leList[0,0]
    y = np.linspace(leList[0,1], leList[1,1], nrib1)
    for k in range(1, nrib1-1):
        e, dist = d.find_closest_edge(x, y[k])
        d.split_face(e.face, x, y[k], droot[0], droot[1])

    x = np.linspace(leList[1,0], leList[2,0], nrib2)
    y = np.linspace(leList[1,1], leList[2,1], nrib2)
    for k in range(1, nrib2-1):
        e, dist = d.find_closest_edge(x[k], y[k])
        d.split_face(e.face, x[k], y[k], drib[0], drib[1])

    # Add all of the points
    x, edge_conn, face_conn, face_sense = d.get_connectivity()
    x = np.array(x, dtype=np.float)
    X = np.zeros((x.shape[0], 2, 3), dtype=np.float)
    X[:,0,:2] = x[:]
    X[:,1,:2] = x[:]
    
    return X, edge_conn, face_conn, face_sense

# Specify the boundary of the planform
# leList = [[25.0, 0.0],
#           [25.0, 3.15],
#           [49.5, 36.0]]
# teList = [[30.3, 0.0],
#           [30.3, 3.15],
#           [50.0, 36.0]]

# Set the leading edge list and the trailing edge lists
leList = [[26.4, 0.0],
          [26.4, 3.0544008],
          [46.19, 29.35]]
teList = [[32.32, 0.0],
          [32.32, 3.0544008],
          [46.87, 29.35]]

# Scale the positions
leList = 25.4*np.array(leList)
teList = 25.4*np.array(teList)

# Create the egads context
ctx = egads.context()

# Set the output file
omlfile = 'ucrm_9_oml.step'

# Load the OML model
oml_body = load_oml_model(ctx, leList, teList, omlfile)
xmin, xmax = oml_body.getBoundingBox()

# Create the ribs/spars
X, edge_conn, face_conn, face_sense = compute_ribspar_edges(leList, teList, 
    nrib1=2, nrib2=4)
X = np.array(X)

# Set the bounds on z so that they encompass the box
X[:,0,2] = xmin[2]-0.1
X[:,1,2] = xmax[2]+0.1

# Create the faces, shell and body
faces = create_faces(ctx, X, edge_conn)

# Create the shell/body object
shell = ctx.makeTopology(egads.SHELL, egads.OPEN, children=faces)
body = ctx.makeTopology(egads.BODY, egads.SHEETBODY, children=[shell])
model = ctx.makeTopology(egads.MODEL, children=[body])

# Save the model
model.saveModel('ribspars.step', overwrite=True)

# Create the CRM step file
create_crm_model(ctx, oml_body, body, filename='ucrm.step')
