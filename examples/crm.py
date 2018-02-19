import numpy as np
from egads4py import egads

crm = egads.pyego()
crm.setOutLevel(2)
crm.loadModel('final_surface.igs')

# Get the 'bodies' from the model
geo, oclass, mtype, lim, bodies, sense = crm.getTopology()

# Create the faces
iges_list = []
for body in bodies:
    faces = body.getBodyTopos(egads.FACE)
    for f in faces:
        iges_list.append(f)

oml = crm.sewFaces(iges_list)

# Loop over the bodies and create a
n1 = np.array([700, 0, -25], dtype=np.float)
n2 = np.array([1200, 700, -25], dtype=np.float)
n3 = np.array([1200, 700, 200], dtype=np.float)
n4 = np.array([700, 0, 200], dtype=np.float)

# Create a box
box = crm.makeSolidBody(egads.BOX, rdata=[n1, [100.0, 100.0, 100.0]])

# # Loop over the bodies and create a
# n1 = np.array([700, 0, -25], dtype=np.float)
# n2 = np.array([1200, 700, -25], dtype=np.float)
# n3 = np.array([1200, 700, 200], dtype=np.float)
# n4 = np.array([700, 0, 200], dtype=np.float)

# oclass = egads.NODE
# node1 = crm.makeTopology(oclass, rdata=n1, refgeo=False)
# node2 = crm.makeTopology(oclass, rdata=n2, refgeo=False)
# node3 = crm.makeTopology(oclass, rdata=n3, refgeo=False)
# node4 = crm.makeTopology(oclass, rdata=n4, refgeo=False)

# # Make the lines between the edges
# oclass = egads.CURVE
# mtype = egads.LINE
# topo_class = egads.EDGE
# topo_type = egads.TWONODE

# d1 = n2 - n1 
# l1 = egads.pyego(crm)
# l1.makeGeometry(oclass, mtype, rdata=[n1, d1])
# e1 = l1.makeTopology(topo_class, topo_type, [node1, node2],
#                      rdata=[0, np.sqrt(np.dot(d1, d1))])

# d2 = n3 - n2 
# l2 = egads.pyego(crm)
# l2.makeGeometry(oclass, mtype, rdata=[n2, d2])
# e2 = l2.makeTopology(topo_class, topo_type, [node2, node3],
#                      rdata=[0, np.sqrt(np.dot(d2, d2))])

# d3 = n4 - n3
# l3 = egads.pyego(crm)
# l3.makeGeometry(oclass, mtype, rdata=[n3, d3])
# e3 = l3.makeTopology(topo_class, topo_type, [node3, node4],
#                      rdata=[0, np.sqrt(np.dot(d3, d3))])

# d4 = n1 - n4
# l4 = egads.pyego(crm)
# l4.makeGeometry(oclass, mtype, rdata=[n4, d4])
# e4 = l4.makeTopology(topo_class, topo_type, [node4, node1],
#                      rdata=[0, np.sqrt(np.dot(d4, d4))])

# # Create a loop for the face
# loop = crm.makeTopology(egads.LOOP, mtype=egads.CLOSED,
#                         children=[e1, e2, e3, e4], sens=[1, 1, 1, 1],
#                         refgeo=False)

# # Create the faces
# oclass = egads.SURFACE
# mtype = egads.PLANE
# plane = egads.pyego(crm)
# plane.makeGeometry(oclass, mtype, rdata=[n1, d1, d2])
# face = plane.makeTopology(egads.FACE, egads.SFORWARD,
#                           children=[loop], sens=[1],
#                           rdata=[0, np.sqrt(np.dot(d1, d1)),
#                                  0, np.sqrt(np.dot(d2, d2))])

inter = box.intersection(oml)

fname = 'test-surface.step'
box.saveModel(fname)
