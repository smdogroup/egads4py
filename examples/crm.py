import numpy as np
import os
import sys
from egads4py import egads
from dcel import dcel

def create_faces(ctx, X, frame_edges, direction):
    # Create the nodes
    nodes = []
    nx = len(X)
    ne = len(frame_edges)

    # Create the bottom nodes
    for i in range(nx):
        oclass = egads.NODE
        node = ctx.makeTopology(oclass, rdata=X[i])
        nodes.append(node)

    # Create the top nodes
    for i in range(nx):
        oclass = egads.NODE
        x = X[i] + direction
        node = ctx.makeTopology(oclass, rdata=x)
        nodes.append(node)

    # Edge list
    edges = []

    # Create the bottom edges
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
        edge.attributeAdd('name', egads.ATTRSTRING, 'bottom')
        edges.append(edge)

    # Create the top edges
    for e in frame_edges:
        x = X[e[0]] + direction
        d = X[e[1]] - X[e[0]] 
        oclass = egads.CURVE
        mtype = egads.LINE
        line = ctx.makeGeometry(oclass, mtype, rdata=[x, d])

        topo_class = egads.EDGE
        topo_type = egads.TWONODE
        edge = ctx.makeTopology(topo_class, topo_type, geom=line,
                                children=[nodes[e[0]+nx], 
                                          nodes[e[1]+nx]],
                                rdata=[0, np.sqrt(np.dot(d, d))])
        edge.attributeAdd('name', egads.ATTRSTRING, 'top')
        edges.append(edge)

    # Create the bottom to top edges
    for i in range(nx):
        d = direction
        oclass = egads.CURVE
        mtype = egads.LINE
        line = ctx.makeGeometry(oclass, mtype, rdata=[X[i], d])

        topo_class = egads.EDGE
        topo_type = egads.TWONODE
        edge = ctx.makeTopology(topo_class, topo_type, geom=line,
                                children=[nodes[i], nodes[nx+i]],
                                rdata=[0, np.sqrt(np.dot(d, d))])
        edges.append(edge)

    # Create the faces
    faces = []
    for i, e in enumerate(frame_edges):
        elist = [edges[i], edges[ne+i], 
                 edges[2*ne+e[0]], edges[2*ne+e[1]]]
        loop, nedge = ctx.makeLoop(elist)
        face = ctx.makeFace(loop, mtype=egads.SFORWARD)
        face.attributeAdd('name', egads.ATTRSTRING, 'box%d'%(i))
        faces.append(face)

    return faces

def create_crm_model(ctx, leList, teList, ribspars, 
                     igesfile, surf_index, filename='ucrm.step'):
    # Load the model 
    crm = ctx.loadModel(igesfile)

    # Get the FACEBODIES from the iges file model and prepare 
    # to sew them together
    body_surfs = crm.getChildren()

    # Keep track of the faces that should be in the final wing-box
    wingbox_faces = []

    # Extract the faces 
    surf_faces = []
    for body in [body_surfs[index] for index in surf_index]:
        for face in body.getBodyTopos(egads.FACE):
            surf_faces.append(face)

    # Create the model by sewing all the faces together
    tol = 1e-5
    for itr in range(10):
        print 'tol = ', tol
        oml = ctx.sewFaces(surf_faces, toler=tol, manifold=False)
        if len(oml.getChildren()) == 1:
            break
        else:
            tol *= np.sqrt(10)
    oml.saveModel('oml.step', overwrite=True)

    # Now 'bodies' should be a single SHEETBODY
    body = oml.getChildren()[0]

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
    tol = 5e-4
    oml = ctx.sewFaces(wingbox_faces, toler=tol, manifold=False)
    oml.saveModel(filename, overwrite=True)
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
    X = np.zeros((x.shape[0], 3), dtype=np.float)
    X[:,:2] = x[:]
    
    return X, edge_conn, face_conn, face_sense

# Specify the boundary of the planform
# leList = [[25.0, 0.0],
#           [25.0, 3.15],
#           [49.5, 36.0]]
# teList = [[30.3, 0.0],
#           [30.3, 3.15],
#           [50.0, 36.0]]

# Set the leading edge list and the trailing edge lists
leList = [[26.4, 1e-10],
          [26.4, 3.0544008],
          [46.19, 29.35]]
teList = [[32.32, 1e-10],
          [32.32, 3.0544008],
          [46.87, 29.35]]

# Scale the positions
leList = 25.4*np.array(leList)
teList = 25.4*np.array(teList)

# Create the ribs/spars
X, edge_conn, face_conn, face_sense = compute_ribspar_edges(leList, teList)
X = np.array(X)

# Set the z-locations
X[:,2] = -50.0
direction = np.array([0, 0, 250.0])

# Create the egads context
ctx = egads.context()

# Create the faces, shell and body
faces = create_faces(ctx, X, edge_conn, direction)
shell = ctx.makeTopology(egads.SHELL, egads.OPEN,
                         children=faces)
body = ctx.makeTopology(egads.BODY, egads.SHEETBODY,
                        children=[shell])

# Save the rib/spar arrangement as a model
model = ctx.makeTopology(egads.MODEL, children=[body])
model.saveModel('ribspars.step', overwrite=True)

# Parameters for the uCRM 13.5
# igesfile = 'ucrm_13_5.iges'
# surf_index = [0, 1, 2, 3, 4, 5]

# Parameter for the CRM
igesfile = 'final_surface.igs'
surf_index = [0, 1, 2, 3, 5, 6]

# Create the CRM step file
create_crm_model(ctx, leList, teList, body, igesfile, surf_index,
                 filename='ucrm.step')
