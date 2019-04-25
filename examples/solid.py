from egads4py import egads

ctx = egads.context()

# Create the bodies that are 
x = [0, 0, 0]
d = [0.25, 0.25, 0.25]
b1 = ctx.makeSolidBody(egads.BOX, rdata=[x, d])

x = [0, 0, 0.25]
d = [0.25, 0.25, 0.75]
b2 = ctx.makeSolidBody(egads.BOX, rdata=[x, d])

x = [0.25, 0, 0]
d = [0.75, 0.25, 0.25]
b3 = ctx.makeSolidBody(egads.BOX, rdata=[x, d])

x = [0, 0.25, 0]
d = [0.25, 0.75, 0.25]
b4 = ctx.makeSolidBody(egads.BOX, rdata=[x, d])

model = b1.solidBoolean(b2, egads.FUSION)
for b in [b3, b4]:
    body = model.getChildren()[0]
    model = body.solidBoolean(b, egads.FUSION)

body = model.getChildren()[0]

# Subtract the cylinders that are the attachment points
x1 = [0.85, 0.0, 0.125]
x2 = [0.85, 0.25, 0.125]
radius = 0.05
c1 = ctx.makeSolidBody(egads.CYLINDER, rdata=[x1, x2, radius])

x1 = [0.125, 0.85, 0.0]
x2 = [0.125, 0.85, 0.25]
radius = 0.05
c2 = ctx.makeSolidBody(egads.CYLINDER, rdata=[x1, x2, radius])

x1 = [0.0, 0.125, 0.85]
x2 = [0.25, 0.125, 0.85]
radius = 0.05
c3 = ctx.makeSolidBody(egads.CYLINDER, rdata=[x1, x2, radius])

for c in [c1, c2, c3]:
    model = body.solidBoolean(c, egads.SUBTRACTION)
    body = model.getChildren()[0]

# Write out to a STEP file
model.saveModel('model.step', overwrite=True)
