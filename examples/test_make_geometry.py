from egads4py import egads

ego = egads.pyego()

# Try to make all of the curves
oclass = egads.CURVE

param = 0.5

# Set the x/y/z locations for the different data to use
loc = [0, 0, 0]
xaxis = [1.0, 0.0, 0.0]
yaxis = [0.0, 1.0, 0.0]
direct = [1, 2, 3]
radius = 4.0
minor = 3.0
major = 9.0

# Make a circle
mtype = egads.CIRCLE
ego.makeGeometry(oclass, mtype, rdata=[loc, xaxis, yaxis, radius])
r, p = ego.getRange()
r = ego.evaluate(param)

# Make a line
mtype = egads.LINE
ego.makeGeometry(oclass, mtype, rdata=[loc, direct])
r, p = ego.getRange()
r = ego.evaluate(param)

# Make an ellipse
mtype = egads.ELLIPSE
ego.makeGeometry(oclass, mtype, rdata=[loc, xaxis, yaxis, major, minor])
r, p = ego.getRange()
r = ego.evaluate(param)

# Make a parabola
mtype = egads.PARABOLA
focus = 3.0
ego.makeGeometry(oclass, mtype, rdata=[loc, xaxis, yaxis, focus])
r, p = ego.getRange()
r = ego.evaluate(param)

# Make a hyperbola
mtype = egads.HYPERBOLA
ego.makeGeometry(oclass, mtype, rdata=[loc, xaxis, yaxis, major, minor])
r, p = ego.getRange()
r = ego.evaluate(param)

# Create a bezier curve
mtype = egads.BEZIER
bitflag = 0 # 2 for rational, 4 for periodic
degree = 3
ncp = 4
pts = [0, 0, 0,  1, 1, 1,  2, 2, 2,  3, 3, 3]
ego.makeGeometry(oclass, mtype, rdata=pts, idata=[bitflag, degree, ncp])
r, p = ego.getRange()
r = ego.evaluate(param)

# Create a bpsline curve
mtype = egads.BSPLINE
bitflag = 0 # 2 for rational, 4 for periodic
degree = 3
ncp = 4
nknots = 4 + degree + 1
knots = [0, 0, 0, 0, 1, 1, 1, 1]
ego.makeGeometry(oclass, mtype, rdata=[knots, pts],
                 idata=[bitflag, degree, ncp, nknots])
r, p = ego.getRange()
ego.evaluate(param)

# Try to make all of the pcurves
oclass = egads.PCURVE

# Reset the x/y axis data
loc = [0, 0]
xaxis = [1.0, 0.0]
yaxis = [0.0, 1.0]
direct = [1, 4]

# Make a circle
mtype = egads.CIRCLE
ego.makeGeometry(oclass, mtype, rdata=[loc, xaxis, yaxis, radius])
ego.evaluate(param)

# Make a line
mtype = egads.LINE
ego.makeGeometry(oclass, mtype, rdata=[loc, direct])
ego.evaluate(param)

# Make an ellipse
mtype = egads.ELLIPSE
ego.makeGeometry(oclass, mtype, rdata=[loc, xaxis, yaxis, major, minor])

# Make a parabola
mtype = egads.PARABOLA
focus = 3.0
ego.makeGeometry(oclass, mtype, rdata=[loc, xaxis, yaxis, focus])
                 
# Make a hyperbola
mtype = egads.HYPERBOLA
ego.makeGeometry(oclass, mtype, rdata=[loc, xaxis, yaxis, major, minor])
ego.evaluate(param)

# Create a bezier curve
mtype = egads.BEZIER
bitflag = 0 # 2 for rational, 4 for periodic
degree = 3
ncp = 4
pts = [0, 0,  1, 1,  2, 2,  3, 3]
ego.makeGeometry(oclass, mtype, rdata=pts, idata=[bitflag, degree, ncp])
ego.evaluate(param)

# Create a bpsline curve
mtype = egads.BSPLINE
bitflag = 0 # 2 for rational, 4 for periodic
degree = 3
ncp = 4
nknots = 4 + degree + 1
knots = [0, 0, 0, 0, 1, 1, 1, 1]
ego.makeGeometry(oclass, mtype, rdata=[knots, pts],
                 idata=[bitflag, degree, ncp, nknots])
ego.evaluate(param)
