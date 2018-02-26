import numpy as np
import os
import sys
from egads4py import egads
from dcel import dcel

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

    body_model = ctx.sewFaces(all_faces, manifold=False)

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

        for face in body.getBodyTopos(egads.FACE):
            attr = face.attributeRet('name')
            if attr == 'box':
                # Check if the edges have
                has_top = True
                has_bottom = True                    
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
                # has_top = False
                # has_bottom = False
                # count = 0
                # for edge in body.getBodyTopos(egads.EDGE, ref=face):
                #     index = body.indexBodyTopo(edge)
                #     if 'box' in edge_attrs[index]:
                #         count += 1
                #     for attr in edge_attrs[index]:
                #         if attr == 'top':
                #             has_top = True
                #         elif attr == 'bottom':
                #             has_bottom = True

                # if not (has_top and has_bottom):
                # wingbox_faces.append(face)
                pass
            
    # Sew the faces together
    print 'Sewing it all together one last time...'
    oml = ctx.sewFaces(wingbox_faces, toler=1e-3, manifold=False)

    oml.saveModel(filename, overwrite=True)
    return

def compute_ribspar_edges():
    # Specify the boundary of the planform
    # leList = [[25.0, 0.0], 
    #           [25.0, 3.15], 
    #           [49.5, 36.0]]
    # teList = [[30.3, 0.0], 
    #           [30.3, 3.15], 
    #           [50.0, 36.0]]

    # # Scale the positions
    # leList = 25.4*np.array(leList)
    # teList = 25.4*np.array(teList)

    leList = [[27.0, 0.0],
              [27.0, 3.15],
              [46.0, 29.0]]
    teList = [[30.5, 0.0],
              [30.5, 3.15],
              [46.5, 29.0]]
    leList = 25.4*np.array(leList)
    teList = 25.4*np.array(teList)

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

    # Set the number of ribs/spars
    nrib1 = 4
    nrib2 = 44

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

# Create the ribs/spars
X, edge_conn, face_conn, face_sense = compute_ribspar_edges()
X = np.array(X)
X[:,2] = -50.0
dist = 250.0
direction = [0, 0, 1]

if 'plot' in sys.argv:
    import matplotlib.pylab as plt
    plt.figure()
    for c in edge_conn:
        plt.plot([X[c[0], 0], X[c[1], 0]], [X[c[0],1], X[c[1],1]])
    plt.axis('equal')
    plt.show()

ctx = egads.context()

# Create all of the edges
edges = create_edges(ctx, X, edge_conn)

edge_loops = []
edge_sense = []
for face_list, sense_list in zip(face_conn, face_sense):
    for e, s in zip(face_list, sense_list):
        if edges[e] is not None:
            edge_loops.append(edges[e])
            edges[e] = None
            edge_sense.append(s)

loop = ctx.makeTopology(egads.LOOP, egads.OPEN,
                        children=edge_loops, sens=edge_sense)
loop, nedges = ctx.makeLoop(edge_loops)

# print nedges
# body = loop.extrude(dist, direction)
# faces = body.getBodyTopos(egads.FACE)

ribspar = ctx.sewFaces(faces, manifold=False)
rsbodies = ribspar.getChildren()

print 'len(rsbodies) = ', len(rsbodies)

# Parameters for the uCRM 13.5
# igesfile = 'ucrm_13_5.iges'
# top_index = [0, 2, 4]
# bottom_index = [1, 3, 5]

# Parameter for the CRM
igesfile = 'final_surface.igs'
top_index = [0, 2, 5]
bottom_index = [1, 3, 6]

# Create the CRM step file
create_crm_model(ctx, rsbodies[0], igesfile, top_index, bottom_index,
                 filename='ucrm.step')
