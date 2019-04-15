from tmr import TMR
import numpy as np

def checkCoincidentVerts(v1, v2, abs_tol=1e-6):
    '''
    Check if two verticies are coincident, within a given tolerance
    '''
    coincident = False

    # Trivial check if they are the same object
    if v1 == v2:
        return True

    # Get each point and compare them
    pt1 = v1.evalPoint()
    pt2 = v2.evalPoint()
    dx = np.abs(pt2-pt1)
    if np.amax(dx) < abs_tol:
        coincident = True

    return coincident

def checkCoincidentEdges(e1, e2, abs_tol=1e-6):
    '''
    Check if two edges are coincident, within a given tolerance
    '''
    coincident = False

    # Trivial check if they are the same object
    if e1 == e2:
        return True

    # Get the verts on each edge
    e1_v1, e1_v2 = e1.getVertices()
    e2_v1, e2_v2 = e2.getVertices()

    # Check if the endpoints are coincident
    if ((checkCoincidentVerts(e1_v1, e2_v1, abs_tol=abs_tol) and
        checkCoincidentVerts(e1_v2, e2_v2, abs_tol=abs_tol)) or
        (checkCoincidentVerts(e1_v2, e2_v1, abs_tol=abs_tol) and
        checkCoincidentVerts(e1_v1, e2_v2, abs_tol=abs_tol))):
        # Check if these are degenerate edges
        if (e1.isDegenerate() and e2.isDegenerate()):
            return True # same and both degenerate
        elif (e1.isDegenerate() + e2.isDegenerate() == 1):
            return False # only one is degenerate
        # Check a midpoint on each
        else:
            e1_tmin, e1_tmax = e1.getRange()
            e1_pt = e1.evalPoint(0.5*(e1_tmin+e1_tmax))
            t2 = e2.invEvalPoint(e1_pt)
            e2_pt = e1.evalPoint(t2)
            dx = np.abs(e2_pt-e1_pt)
            if np.amax(dx) < abs_tol:
                coincident = True

    return coincident

def checkCoincidentEdgeLoops(eloop1, eloop2, abs_tol=1e-6):
    '''
    Check if two edge loops are coincident, within a given tolerance
    '''
    coincident = False

    # Trivial check if they are the same object
    if eloop1 == eloop2:
        return True

    edges1, dirs1 = eloop1.getEdgeLoop()
    edges2, dirs2 = eloop2.getEdgeLoop()

    # If they have different number of edges, not coincident
    if len(edges1) != len(edges2):
        return False

    # Compare edge loops, comparing both directions
    nmatches = 0
    for e1 in edges1:
        for e2 in edges2:
            if checkCoincidentEdges(e1, e2, abs_tol=abs_tol):
                nmatches += 1

    return (nmatches == len(edges1))

def checkCoincidentFaces(f1, f2, abs_tol=1e-6):
    '''
    Check if two faces are coincident, within a given tolerance
    '''
    coincident = False

    # Trivial check if they are the same object
    if f1 == f2:
        return True

    nloop1 = f1.getNumEdgeLoops()
    nloop2 = f2.getNumEdgeLoops()
    # No match if the faces have different numbers of edge loops
    if nloop1 != nloop2:
        return False

    # Check that each edge loop has a match
    nmatched = 0
    for i in range(nloop1):
        l1 = f1.getEdgeLoop(i)
        for j in range(nloop2):
            l2 = f2.getEdgeLoop(j)
            if checkCoincidentEdgeLoops(l1, l2, abs_tol=abs_tol):
                nmatched += 1
                break

    if nmatched == nloop1:
        coincident = True

    return coincident

def findMatchingFaces(geo1, geo2):
    '''
    Take in two geometries, find the matching faces, and set them as copies
    '''
    faces1 = geo1.getFaces()
    faces2 = geo2.getFaces()
    for f1 in faces1:
        for f2 in faces2:
            if checkCoincidentFaces(f1, f2):
                setCopyFaces(f1, f2)

    return

def setCopyFaces(source, target, abs_tol=1e-6):
    '''
    Take in two coincident faces, set the target as a copy of the source,
    find their coincident edges, and set the edges as copies as well
    '''
    target.setCopySource(source, orient=1)
    
    source_list = []
    for i in range(source.getNumEdgeLoops()):
        sloop = source.getEdgeLoop(i)     
        s_edges, dirs = sloop.getEdgeLoop()
        source_list.extend(s_edges)

    target_list = []
    for j in range(target.getNumEdgeLoops()):
        tloop = target.getEdgeLoop(j)
        t_edges, dirs = tloop.getEdgeLoop()
        target_list.extend(t_edges)

    for s in source_list:
        for t in target_list:
            if checkCoincidentEdges(t, s, abs_tol=abs_tol):
                # Remove the target from the list of targets
                target_list.remove(t)
                sv1, sv2 = s.getVertices()
                tv1, tv2 = t.getVertices()

                # Set the target edge source edge
                t.setCopySource(s)
                if checkCoincidentVerts(sv1, tv1, abs_tol=abs_tol):
                    tv1.setCopySource(sv1)
                    tv2.setCopySource(sv2)
                else:
                    tv1.setCopySource(sv2)
                    tv2.setCopySource(sv1)
                break
            
    return
    
