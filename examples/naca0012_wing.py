import numpy as np
from egads4py import egads

def dist(y1, y2):
    return np.sqrt(np.dot(y2 - y1, y2 - y1))

ctx = egads.context()
ctx.setOutLevel(0)

# Set the overall parameters of the wing, including 1/4 chord sweep,
# taper ratio, aspect ratio etc.
chord_sweep = 15.0*np.pi/180.0
root_chord = 4.0
taper_ratio = 0.35
AR = 4.0
root_le_fraction = 0.125
root_te_fraction = 0.75
tip_le_fraction = 0.15
tip_te_fraction = 0.55

# Compute additional geometric parameters
tip_chord = taper_ratio*root_chord
b = AR*root_chord*(1.0 + taper_ratio)
semi_span = b/2.0
S = b**2/AR

# Load the foil data
Xfoil = np.loadtxt('naca0012.dat')
xlen = Xfoil.shape[0]
xlen = xlen//2
Xorig = Xfoil[:xlen,:]
Xfoil[:xlen,:] = Xorig[::-1,:]

# Set the locations of the non-dimensional airfoil
X = np.zeros((Xfoil.shape[0], 3))
X[:,0] = Xfoil[:,0]
X[:,2] = Xfoil[:,1]

# Make the trailing edge sharp
X[:,2] += np.linspace(-X[0,2], -X[-1,2], X.shape[0])

# Set airfoil coordinates at the root
Xroot = root_chord*X
Xroot[:,0] -= 0.25*root_chord

# Set the airfoil coordinates at the tip
Xtip = tip_chord*X
Xtip[:,0] -= 0.25*tip_chord
Xtip[:,0] += semi_span*np.tan(chord_sweep)
Xtip[:,1] = semi_span

# Find the root/tip foil coordinates via a b-spline fit in egads
root_curve = ctx.approximate(Xroot, maxdeg=0, tol=1e-6)
root_range, periodic = root_curve.getRange()
Xroot_node = Xroot[0,:]
root_node = ctx.makeTopology(egads.NODE, rdata=Xroot_node)
root_edge = ctx.makeTopology(egads.EDGE, mtype=egads.ONENODE,
                             geom=root_curve, children=[root_node],
                             rdata=root_range)
root_loop = ctx.makeTopology(egads.LOOP, mtype=egads.CLOSED,
                             children=[root_edge], sens=[1])
root_face = ctx.makeFace(root_loop, egads.SREVERSE)

tip_curve = ctx.approximate(Xtip, maxdeg=0, tol=1e-6)
tip_range, periodic = tip_curve.getRange()
Xtip_node = Xtip[0,:]
tip_node = ctx.makeTopology(egads.NODE, rdata=Xtip_node)
tip_edge = ctx.makeTopology(egads.EDGE, mtype=egads.ONENODE,
                            geom=tip_curve, children=[tip_node],
                            rdata=tip_range)
tip_loop = ctx.makeTopology(egads.LOOP, mtype=egads.CLOSED,
                            children=[tip_edge], sens=[1])
tip_face = ctx.makeFace(tip_loop, egads.SREVERSE)

# Create the semi-span wing
semi_wing = ctx.ruled([tip_face, root_face])
semi_wing = ctx.makeTopology(egads.MODEL, children=[semi_wing])

# Create the block that we will extrude through the wing
oclass = egads.NODE
zloc = -root_chord
x1 = [root_chord*(root_le_fraction - 0.25), 0, zloc]
x2 = [root_chord*(root_te_fraction - 0.25), 0, zloc]
x3 = [tip_chord*(tip_le_fraction - 0.25) + semi_span*np.tan(chord_sweep),
      semi_span, zloc]
x4 = [tip_chord*(tip_te_fraction - 0.25) + semi_span*np.tan(chord_sweep),
      semi_span, zloc]

x1 = np.array(x1)
x2 = np.array(x2)
x3 = np.array(x3)
x4 = np.array(x4)

n1 = ctx.makeTopology(oclass, rdata=x1)
n2 = ctx.makeTopology(oclass, rdata=x2)
n3 = ctx.makeTopology(oclass, rdata=x3)
n4 = ctx.makeTopology(oclass, rdata=x4)

oclass = egads.CURVE
mtype = egads.LINE
topo_class = egads.EDGE
topo_type = egads.TWONODE
l1 = ctx.makeGeometry(oclass, mtype, rdata=[x1, x2 - x1])
e1 = ctx.makeTopology(topo_class, topo_type, geom=l1,
                      children=[n1, n2],
                      rdata=[0, dist(x1, x2)])
l2 = ctx.makeGeometry(oclass, mtype, rdata=[x3, x4 - x3])
e2 = ctx.makeTopology(topo_class, topo_type, geom=l2,
                      children=[n3, n4],
                      rdata=[0, dist(x3, x4)])
l3 = ctx.makeGeometry(oclass, mtype, rdata=[x1, x3 - x1])
e3 = ctx.makeTopology(topo_class, topo_type, geom=l3,
                      children=[n1, n3],
                      rdata=[0, dist(x1, x3)])
l4 = ctx.makeGeometry(oclass, mtype, rdata=[x2, x4 - x2])
e4 = ctx.makeTopology(topo_class, topo_type, geom=l4,
                      children=[n2, n4],
                      rdata=[0, dist(x2, x4)])

# Create the loop
elist = [e1, e4, e2, e3]
senses = [1, 1, -1, -1]
loop = ctx.makeTopology(egads.LOOP, egads.CLOSED,
                        children=elist, sens=senses)

# Create the face from the loop and a body from an extrusion
face = ctx.makeFace(loop, mtype=egads.SFORWARD)
body = face.extrude(2*root_chord, [0, 0, 1])

semi_wing = semi_wing.solidBoolean(body, egads.INTERSECTION)
body = semi_wing.getChildren()[0]
shell = body.getChildren()[0]
faces = shell.getChildren()
faces[0].attributeAdd('name', egads.ATTRSTRING, 'bottom')
faces[1].attributeAdd('name', egads.ATTRSTRING, 'tip')
faces[2].attributeAdd('name', egads.ATTRSTRING, 'leading_spar')
faces[3].attributeAdd('name', egads.ATTRSTRING, 'root')
faces[4].attributeAdd('name', egads.ATTRSTRING, 'trailing_spar')
faces[5].attributeAdd('name', egads.ATTRSTRING, 'top')

semi_wing.saveModel('semi_wing.step', overwrite=True)
semi_wing.saveModel('semi_wing.egads', overwrite=True)