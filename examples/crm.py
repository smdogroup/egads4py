import numpy as np
import os
import sys
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

def create_crm_model(ctx, ribspars, igesfile, top_index, bottom_index,
                     filename='ucrm.step'):
    # Load the model 
    crm = ctx.loadModel(igesfile)

    # Get the FACEBODES from the iges file model and prepare 
    # to sew them together
    geo, oclass, mtype, lim, body_surfs, sense = crm.getTopology()

    top_surfs = [body_surfs[index] for index in top_index]
    bottom_surfs = [body_surfs[index] for index in bottom_index]

    # Extract the faces and sew them
    all_faces = []
    for body in top_surfs:
        faces = body.getBodyTopos(egads.FACE)
        for face in faces:
            face.attributeAdd('name', egads.ATTRSTRING, 'top') 
            all_faces.append(face)    

    for body in bottom_surfs:
        faces = body.getBodyTopos(egads.FACE)
        for face in faces:
            face.attributeAdd('name', egads.ATTRSTRING, 'bottom') 
            all_faces.append(face)    

    print 'Sewing iges surfaces together into a single body...'
    body_model = ctx.sewFaces(all_faces, toler=1e-3, manifold=False)

    # Now 'bodies' should be a single SHEETBODY
    bodies = body_model.getChildren()

    for face in ribspars.getBodyTopos(egads.FACE):
        face.attributeAdd('name', egads.ATTRSTRING, 'box')

    new_bodies = []
    for body in bodies:
        print 'Intersecting iges surface body...'
        wirebody, pairs = body.intersection(ribspars)
        if wirebody is not None:
            print 'Impriting iges surface body...'
            new_body = body.imprintBody(pairs)
            new_bodies.append(new_body)
        
        print 'Intersecting ribspar body...'
        wirebody, pairs = ribspars.intersection(body)
        if wirebody is not None:
            ribspars = ribspars.imprintBody(pairs)

    # Now determine which new surfaces to keep and which to discard
    new_bodies.append(ribspars)
    all_faces = []
    for body in new_bodies:
        all_faces.extend(body.getBodyTopos(egads.FACE))

    # Sew all the faces together
    print 'Sewing all the surfaces together into a single body...'
    everything = ctx.sewFaces(all_faces, toler=1e-3, manifold=False)

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

        for face in body.getBodyTopos(egads.FACE):
            attr = face.attributeRet('name')
            if attr == 'box':
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
                    wingbox_faces.append(face)
            elif attr == 'top' or attr == 'bottom':
                has_top = False
                has_bottom = False
                count = 0
                for edge in body.getBodyTopos(egads.EDGE, ref=face):
                    index = body.indexBodyTopo(edge)
                    if 'box' in edge_attrs[index]:
                        count += 1
                    for attr in edge_attrs[index]:
                        if attr == 'top':
                            has_top = True
                        elif attr == 'bottom':
                            has_bottom = True

                if not (has_top and has_bottom):
                    wingbox_faces.append(face)
            
    # Sew the faces together
    print 'Sewing it all together one last time...'
    oml = ctx.sewFaces(wingbox_faces, toler=1e-3, manifold=False)

    oml.saveModel(filename, overwrite=True)
    return

def compute_ribspar_edges():
    # Specify the boundary of the planform
    leList = [[25.0, 0.0, 0.0], 
              [25.0, 3.15, 0.0], 
              [49.5, 35.95, 0.0]]
    teList = [[30.3, 0.0, 0.0], 
              [30.3, 3.15, 0.0], 
              [50.0, 35.95, 0.0]]

    # Scale the Xpts
    leList = np.array(leList)
    teList = np.array(teList)
    leList *= 25.4
    teList *= 25.4

    # Set the number of ribs/spars
    nrib1 = 3
    nrib2 = 45
    nribs = nrib1 + nrib2

    # Set the leading edge points
    le_ypts = np.linspace(leList[0,1], leList[1,1], nrib1+1)[:-1]
    le_ypts = np.append(le_ypts, 
                        np.linspace(leList[1,1], leList[2,1], nrib2))

    le_xpts = np.linspace(leList[0,0], leList[1,0], nrib1+1)[:-1]
    le_xpts = np.append(le_xpts,
                        np.linspace(leList[1,0], leList[2,0], nrib2))

    # Set the trailing edge points
    te_ypts = np.linspace(teList[0,1], teList[1,1], nrib1+1)[:-1]
    te_ypts = np.append(te_ypts, 
                        np.linspace(teList[1,1], teList[2,1], nrib2))

    te_xpts = np.linspace(teList[0,0], teList[1,0], nrib1+1)[:-1]
    te_xpts = np.append(te_xpts,
                        np.linspace(teList[1,0], teList[2,0], nrib2))

    # Add all of the points
    X = []
    for k in range(nribs):
        X.append([le_xpts[k], le_ypts[k], 0.0])
    for k in range(nribs):
        X.append([te_xpts[k], te_ypts[k], 0.0])

    conn = []
    for k in range(nribs-1):
        conn.append([k+1, k+1 + nribs])
    for k in range(nribs-1):
        conn.append([k, k+1])
        conn.append([k+nribs, k+1+nribs])
        
    return X, conn 

# Create the ribs/spars
X, conn = compute_ribspar_edges()
X = np.array(X)
X[:,2] = -50.0
dist = 250.0
direction = [0, 0, 1]

if 'plot' in sys.argv:
    import matplotlib.pylab as plt
    plt.figure()
    for c in conn:
        plt.plot([X[c[0], 0], X[c[1], 0]], [X[c[0],1], X[c[1],1]])
    plt.axis('equal')
    plt.show()

ctx = egads.context()

# Create all of the edges
edges = create_edges(ctx, X, conn)

# Create the edges that form the ribs/spars
nedges = len(edges)
faces = []
while nedges > 0:
    loop, nedges = ctx.makeLoop(edges)
    body = loop.extrude(dist, direction)
    faces.extend(body.getBodyTopos(egads.FACE))

ribspar = ctx.sewFaces(faces, toler=1e-5, manifold=False)
ribspar.saveModel('ribspar.step', overwrite=True)
rsbodies = ribspar.getChildren()

# Set the source file and the top/bottom surfaces
igesfile = 'ucrm_13_5.iges'
top_index = [0, 2, 4]
bottom_index = [1, 3, 5]

# Create the CRM step file
create_crm_model(ctx, rsbodies[0], igesfile, top_index, bottom_index,
                 filename='ucrm.step')
